import json
import os
from pathlib import Path
from unittest.mock import patch

import pytest

import qcelemental as qcel

from qcengine.testing import _qcschema_data_path


def _example_files():
    return sorted(_qcschema_data_path.glob("v*/*/*.json"))


files = _example_files()
params = files or [pytest.param(None, id="no-generated-examples")]


def _model_and_version(path: Path):
    version_name, model_name = path.relative_to(_qcschema_data_path).parts[:2]
    version = int(version_name.removeprefix("v"))
    namespace = {1: qcel.models.v1, 2: qcel.models.v2}[version]
    return getattr(namespace, model_name), version


def _schema(model, version):
    if version == 1:
        return model.schema()

    try:
        return model.model_json_schema()
    except KeyError as exc:
        if exc.args != ("dtype",):
            raise

        # QCElemental 0.50.4 cannot find array dtype metadata after Pydantic
        # 2.13 wraps its core schema. Use QCElemental's own schema-export
        # fallback until the recursive dtype lookup is available in a release.
        with patch.dict(os.environ, {"SPHINX_BUILD": "1"}):
            return model.model_json_schema()


def _validate_json_schema(instance, model, version):
    import jsonschema

    schema = _schema(model, version)
    validator_class = jsonschema.validators.validator_for(schema)
    validator_class.check_schema(schema)
    validator_class(schema).validate(instance)


@pytest.mark.parametrize("fl", params, ids=lambda fl: str(fl.relative_to(_qcschema_data_path)) if fl else None)
def test_qcschema_example(fl, request):
    if not request.config.getoption("--validate-qcschema-examples", default=False):
        pytest.skip("QCSchema examples are checked only with --validate-qcschema-examples")
    if fl is None:
        pytest.fail("No generated QCSchema examples found; run pytest --qcschema-examples first")

    model, version = _model_and_version(fl)
    raw = fl.read_text()
    data = json.loads(raw)

    instance = model.parse_raw(raw) if version == 1 else model.model_validate_json(raw)
    _validate_json_schema(data, model, version)

    target_version = {1: 2, 2: 1}[version]
    converted = instance.convert_v(target_version)
    converted_model = type(converted)
    converted_data = json.loads(converted.model_dump_json(exclude_unset=True, exclude_none=True))
    _validate_json_schema(converted_data, converted_model, target_version)
