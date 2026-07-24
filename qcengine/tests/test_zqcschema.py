import json
from functools import cache
from pathlib import Path

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


def _find_array_dtype(schema):
    if isinstance(schema, dict):
        metadata = schema.get("metadata", {})
        if "dtype" in metadata:
            return metadata["dtype"]
        for value in schema.values():
            dtype = _find_array_dtype(value)
            if dtype is not None:
                return dtype
    elif isinstance(schema, (list, tuple)):
        for value in schema:
            dtype = _find_array_dtype(value)
            if dtype is not None:
                return dtype
    return None


def _restore_array_dtype_metadata(schema, restored):
    if isinstance(schema, dict):
        metadata = schema.get("metadata", {})
        hooks = (
            *metadata.get("pydantic_js_functions", ()),
            *metadata.get("pydantic_js_annotation_functions", ()),
        )
        is_array_hook = any(
            getattr(getattr(hook, "__self__", None), "__name__", None) == "ValidatableArrayAnnotation"
            for hook in hooks
        )
        if is_array_hook and "dtype" not in metadata:
            dtype = _find_array_dtype(schema)
            if dtype is not None:
                metadata["dtype"] = dtype
                restored.append(metadata)
        for value in schema.values():
            _restore_array_dtype_metadata(value, restored)
    elif isinstance(schema, (list, tuple)):
        for value in schema:
            _restore_array_dtype_metadata(value, restored)


@cache
def _schema(model, version):
    if version == 1:
        return model.schema()

    try:
        return model.model_json_schema()
    except KeyError as exc:
        if exc.args != ("dtype",):
            raise

    # QCElemental 0.50.4's array JSON-Schema hook expects dtype metadata
    # on a wrapper where Pydantic 2.13 no longer preserves it. Copy each
    # nested dtype to that wrapper for schema export, then restore the model.
    restored = []
    _restore_array_dtype_metadata(model.__pydantic_core_schema__, restored)
    try:
        return model.model_json_schema()
    finally:
        for metadata in restored:
            metadata.pop("dtype")


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
