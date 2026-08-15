(async function () {
  function sleep(ms) {
    return new Promise((r) => setTimeout(r, ms));
  }

  // Wait for RTD sidebar to exist (theme may build it after our script loads)
  async function waitForSidebar(maxTries = 50, delayMs = 100) {
    for (let i = 0; i < maxTries; i++) {
      const search = document.querySelector(".wy-side-nav-search");
      const scroll = document.querySelector(".wy-side-scroll");
      if (search || scroll) return { search, scroll };
      await sleep(delayMs);
    }
    return { search: null, scroll: null };
  }

  const { search, scroll } = await waitForSidebar();
  if (!search && !scroll) {
    console.warn("[version-switcher] RTD sidebar not found.");
    return;
  }

  // Avoid duplicates (important if theme re-renders and script runs again)
  if (document.querySelector(".version-switcher")) return;

  // --- UI ---
  const container = document.createElement("div");
  container.className = "version-switcher";
  container.style.padding = "0.5rem 1rem";
  container.style.display = "flex";
  container.style.gap = "0.5rem";
  container.style.alignItems = "center";
  container.style.justifyContent = "center";

  const label = document.createElement("span");
  label.textContent = "Version:";
  label.style.fontSize = "0.9em";
  label.className = "vs-label";

  const select = document.createElement("select");
  select.setAttribute("aria-label", "Select documentation version");
  select.style.width = "50%";
  select.style.borderRadius = "20px";

  container.appendChild(label);
  container.appendChild(select);

  // Insert under the search box if possible
  if (search && search.parentNode) {
    search.insertAdjacentElement("afterend", container);
  } else {
    scroll.prepend(container);
  }

  // --- Load versions.json ---
  select.disabled = true;
  select.innerHTML = `<option>Loading…</option>`;

  async function fetchJson(url) {
    const resp = await fetch(url, { cache: "no-store" });
    if (!resp.ok) throw new Error(`${resp.status} ${resp.statusText}`);
    return resp.json();
  }

  // For local http://localhost:8000/ this will work.
  // For GitHub project pages it will also work if versions.json is at the root of the published site.
  let versions;
  try {
    versions = await fetchJson(`${window.location.origin}/versions.json`);
  } catch (e) {
    // GitHub project page fallback: /QCElemental/versions.json
    const parts = window.location.pathname.split("/").filter(Boolean);
    if (parts.length) {
      try {
        versions = await fetchJson(`${window.location.origin}/${parts[0]}/versions.json`);
      } catch (e2) {
        console.warn("[version-switcher] Could not load versions.json", e, e2);
        select.innerHTML = `<option>No versions.json</option>`;
        return;
      }
    } else {
      console.warn("[version-switcher] Could not load versions.json", e);
      select.innerHTML = `<option>No versions.json</option>`;
      return;
    }
  }

  // Populate
  select.innerHTML = "";
  for (const v of versions) {
    const opt = document.createElement("option");
    opt.value = v.path;        // "dev/" etc.
    opt.textContent = v.label; // "dev"
    select.appendChild(opt);
  }

  // Determine current "version folder" (dev / v0.30.1 / etc.)
  // If you're at /dev/index.html => currentVersion = "dev"
  // If you're at /index.html => currentVersion = ""
  const pathParts = window.location.pathname.split("/").filter(Boolean);
  const currentVersion = pathParts.length >= 2 ? pathParts[1] : pathParts[0] || "";

  for (const opt of select.options) {
    if (opt.value.replace(/\/+$/, "") === currentVersion) {
      opt.selected = true;
      break;
    }
  }

  select.disabled = false;

  // Navigate on change; for local dev it’s simplest to go to version root
  select.addEventListener("change", () => {
    const target = select.value; // e.g. "dev/"
    // If on GitHub project site, keep the project prefix automatically.
    const parts = window.location.pathname.split("/").filter(Boolean);
    const prefix = parts.length ? `/${parts[0]}/` : "/";
    window.location.href = `${window.location.origin}${prefix}${target}`;
  });
})();

