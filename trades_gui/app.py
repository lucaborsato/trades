"""Web GUI for creating TRADES configuration directories."""
import os
import sys

import streamlit as st

sys.path.insert(0, os.path.dirname(__file__))
from config_generator import (  # noqa: E402
    ANALYSIS_PARAMS,
    ARG_PARAMETERS,
    DEFAULT_CONFIG,
    EMCEE_PARAMS,
    FIT_PARAMETERS,
    NB_COMMENT,
    OC_PARAMS,
    PYDE_PARAMS,
    RUN_PARAMS,
    RV_PARAMS,
    generate_arg_in,
    generate_bodies_lst,
    generate_configuration_yml,
    generate_lm_opt,
    generate_pikaia_opt,
    generate_planet_dat,
    generate_priors,
    generate_pso_opt,
    generate_rv_observations,
    generate_star_dat,
    generate_transit_observations,
    get_planet_label,
    load_example,
    write_system,
)


st.set_page_config(page_title="TRADES Config Generator", layout="wide")
st.markdown("""
<style>
.stButton > button { width: 100%; }
textarea, code { font-family: monospace; }
</style>
""", unsafe_allow_html=True)


def float_input(label, value, key, min_value=None):
    """Free-form numeric entry: no spinner buttons and no rounding of the typed value."""
    if key not in st.session_state:
        st.session_state[key] = repr(float(value))
    text = st.text_input(label, key=key, label_visibility="collapsed")
    cleaned = text.strip()
    if "," in cleaned and "." not in cleaned:
        cleaned = cleaned.replace(",", ".")
    try:
        number = float(cleaned)
    except ValueError:
        number = float(value)
    if min_value is not None and number < min_value:
        number = float(min_value)
    return number


INPUT_KEY_TOKENS = (
    "star_", "_value_", "_min_", "_max_", "_prior_value_", "_sigma_",
    "tt_", "ett_", "dur_", "edur_", "ep_", "lc_", "rv_jd_", "rv_", "erv_",
    "duration_", "nobs_", "rv_count", "rvset_", "label_", "transit_",
    "planet_count", "sys_name", "export_dir", "arg_", "cfg_",
)


def clear_input_state():
    """Drop widget state so inputs re-seed from the current system data."""
    for key in list(st.session_state.keys()):
        if key.startswith("$") or key == "system":
            continue
        if any(token in key for token in INPUT_KEY_TOKENS):
            del st.session_state[key]


@st.dialog("File preview")
def show_preview(filename, content):
    st.caption(filename)
    st.code(content, language=None)


def new_planet(index):
    """Return defaults for one planet."""
    values = {
        "label": get_planet_label(index), "transiting": True,
        "mass": 0.01, "mass_min": 0.0, "mass_max": 20.0,
        "radius": 0.1, "radius_min": 0.0, "radius_max": 20.0,
        "period": 10.0, "period_min": 1.0, "period_max": 1000.0,
        "sma": 999.0, "sma_min": 0.0, "sma_max": 0.0,
        "ecc": 0.0, "ecc_min": 0.0, "ecc_max": 1.0,
        "omega": 90.0, "omega_min": 0.0, "omega_max": 360.0,
        "M0": 0.0, "M0_min": 0.0, "M0_max": 360.0,
        "tp": 9e8, "tp_min": 0.0, "tp_max": 0.0,
        "inc": 90.0, "inc_min": 0.0, "inc_max": 180.0,
        "lnode": 180.0, "lnode_min": 0.0, "lnode_max": 360.0,
        "transit_obs": [], "include_duration": True,
    }
    for key, _, _ in FIT_PARAMETERS:
        values[f"fit_{key}"] = True
        values[f"prior_{key}"] = False
        values[f"prior_{key}_value"] = values[key]
        values[f"prior_{key}_sigma"] = 0.0
    return values


def default_system():
    return {
        "name": "", "n_planets": 2,
        "star": {"mass": 1.0, "mass_err": 0.01, "radius": 1.0,
                 "radius_err": 0.01, "fit_mass": False, "fit_radius": False,
                 "mass_prior": False, "radius_prior": False},
        "planets": [new_planet(0), new_planet(1)],
        "rv_obs": [], "export_dir": "",
        "arg_in": {}, "config": {},
    }


def ensure_planets(system):
    while len(system["planets"]) < system["n_planets"]:
        system["planets"].append(new_planet(len(system["planets"])))
    system["planets"] = system["planets"][:system["n_planets"]]
    for i, planet in enumerate(system["planets"]):
        planet.setdefault("label", get_planet_label(i))
        for key, _, _ in FIT_PARAMETERS:
            planet.setdefault(f"fit_{key}", True)
            planet.setdefault(f"prior_{key}", False)
            planet.setdefault(f"prior_{key}_value", planet.get(key, 0.0))
            planet.setdefault(f"prior_{key}_sigma", 0.0)


def ensure_arg_defaults(args):
    for key, _ptype, default, _comment, _choices in ARG_PARAMETERS:
        args.setdefault(key, default)


def ensure_config_defaults(cfg):
    for full_key, default in DEFAULT_CONFIG.items():
        cfg.setdefault(full_key, default)


def render_params(cfg, prefix, params):
    """Render one configuration section as Parameter | Value | Description rows."""
    header = st.columns([1.5, 1.5, 5])
    for col, text in zip(header, ["Parameter", "Value", "Description"]):
        col.markdown(f"**{text}**")
    for key, ptype, default, comment, choices in params:
        skey = f"cfg_{prefix}_{key}"
        row = st.columns([1.5, 1.5, 5])
        row[0].write(key)
        current = cfg.get(skey, default)
        with row[1]:
            if ptype == "bool":
                value = st.checkbox("", value=bool(current), key=skey, label_visibility="collapsed")
            elif ptype == "choice":
                index = choices.index(current) if current in choices else 0
                value = st.selectbox("", choices, index, key=skey, label_visibility="collapsed")
            elif ptype == "int":
                value = int(st.number_input("", value=int(current), key=skey, label_visibility="collapsed"))
            elif ptype in ("float", "time"):
                value = float_input("", float(current), skey)
            else:
                value = st.text_input("", value="" if current is None else str(current),
                                      key=skey, label_visibility="collapsed")
        row[2].caption(comment or "")
        cfg[skey] = value


def build_file_contents(system):
    """Build the content of every output file from the current system state."""
    star = system["star"]
    planets = system["planets"]
    files = {
        "bodies.lst": generate_bodies_lst(system["n_planets"], star.get("fit_mass", False),
                                          star.get("fit_radius", False), planets),
        "star.dat": generate_star_dat(star.get("mass", 1.0), star.get("mass_err", 0.01),
                                      star.get("radius", 1.0), star.get("radius_err", 0.01), None),
        "arg.in": generate_arg_in(system.get("arg_in", {})),
        "priors.in": generate_priors(planets, star),
        "configuration.yml": generate_configuration_yml(system),
        "lm.opt": generate_lm_opt({}),
        "pikaia.opt": generate_pikaia_opt({}),
        "pso.opt": generate_pso_opt({}),
    }
    for i, planet in enumerate(planets):
        files[f"{planet.get('label', get_planet_label(i))}.dat"] = generate_planet_dat(planet)
    for i, planet in enumerate(planets):
        if planet.get("transiting") and planet.get("transit_obs"):
            files[f"NB{i + 2}_observations.dat"] = generate_transit_observations(
                planet["transit_obs"], planet.get("include_duration", True))
    if system.get("rv_obs"):
        files["obsRV.dat"] = generate_rv_observations(system["rv_obs"])
    return files


def edit_transit_times(planet, body_id):
    """Render and store one planet's NB observations."""
    if not planet.get("transiting", True):
        st.info(f"Enable Transiting to create NB{body_id}_observations.dat.")
        return
    st.markdown(f"**NB{body_id}_observations.dat**")
    st.caption("This interface can be used to create an NB{}_observations.dat file, but you should populate and verify the file yourself before using it with TRADES.".format(body_id))
    planet["include_duration"] = st.checkbox(
        "Include Duration / err_Duration columns",
        value=planet.get("include_duration", True), key=f"duration_{body_id}",
    )
    observations = planet.get("transit_obs", [])
    count = st.number_input("Number of observations", 0, 500, len(observations), key=f"nobs_{body_id}")
    headers = ["Epoch", "Transit_time(days)", "err_Transit_time(days)"]
    if planet["include_duration"]:
        headers.extend(["Duration(min)", "err_Duration(min)"])
    headers.append("LC_id")
    header_cols = st.columns(len(headers))
    for column, header in zip(header_cols, headers):
        column.markdown(f"**{header}**")
    updated = []
    for index in range(count):
        old = observations[index] if index < len(observations) else {}
        cols = st.columns(len(headers))
        epoch = cols[0].number_input("Epoch", value=int(old.get("epoch", index)), key=f"ep_{body_id}_{index}", label_visibility="collapsed")
        with cols[1]:
            tt = float_input("", old.get("tt", 2455000.0), f"tt_{body_id}_{index}")
        with cols[2]:
            ett = float_input("", old.get("ett", 0.001), f"ett_{body_id}_{index}")
        offset = 3
        dur = old.get("dur", 2.5)
        edur = old.get("edur", 0.05)
        if planet["include_duration"]:
            with cols[offset]:
                dur = float_input("", dur, f"dur_{body_id}_{index}")
            with cols[offset + 1]:
                edur = float_input("", edur, f"edur_{body_id}_{index}")
            offset += 2
        lc_id = cols[offset].number_input("LC id", value=int(old.get("lc_id", 1)), key=f"lc_{body_id}_{index}", label_visibility="collapsed")
        updated.append({"epoch": epoch, "tt": tt, "ett": ett, "dur": dur, "edur": edur, "lc_id": lc_id})
    planet["transit_obs"] = updated


if "system" not in st.session_state:
    st.session_state.system = default_system()
system = st.session_state.system

st.title("TRADES Config Generator")
st.caption("Create a complete TRADES input directory in the browser.")

with st.sidebar:
    st.header("Navigation")
    page = st.radio("Section", ["Overview", "Star", "Planets", "RV Observations",
                                "Integration Settings", "Analysis Settings"],
                    label_visibility="collapsed")
    st.divider()
    system["export_dir"] = st.text_input("Export directory", value=system.get("export_dir", ""), key="export_dir")
    if st.button("New System"):
        st.session_state.system = default_system()
        clear_input_state()
        st.rerun()
    example_dir = st.text_input("Load example directory", value="")
    if st.button("Load Example"):
        if os.path.isdir(example_dir):
            st.session_state.system = load_example(example_dir)
            ensure_planets(st.session_state.system)
            clear_input_state()
            st.rerun()
        else:
            st.error("Example directory was not found.")
    st.divider()
    st.subheader("Export")
    if st.button("Generate TRADES files", type="primary"):
        destination = system.get("export_dir", "")
        name = system.get("name", "my_system") or "my_system"
        if not destination:
            st.error("Set an export directory in the sidebar first.")
        else:
            output = os.path.join(destination, name)
            created = write_system(output, system)
            st.success(f"Created {len(created)} files in {output}")
            st.code("\n".join(sorted(created)))

if page == "Overview":
    st.subheader("System Overview")
    system["name"] = st.text_input("System name", value=system.get("name", ""), key="sys_name")
    count = st.number_input("Number of planets", 1, 26, system["n_planets"], key="planet_count")
    if count != system["n_planets"]:
        system["n_planets"] = count
        ensure_planets(system)
        st.rerun()
    st.info(f"The output contains one star and {count} planet(s).")

    ensure_planets(system)
    files = build_file_contents(system)
    order = ["bodies.lst", "star.dat"]
    for i, planet in enumerate(system["planets"]):
        order.append(f"{planet.get('label', get_planet_label(i))}.dat")
    for i, planet in enumerate(system["planets"]):
        if planet.get("transiting") and planet.get("transit_obs"):
            order.append(f"NB{i + 2}_observations.dat")
    if system.get("rv_obs"):
        order.append("obsRV.dat")
    order += ["arg.in", "priors.in", "configuration.yml", "lm.opt", "pikaia.opt", "pso.opt"]

    st.caption("Live preview of the files that Export will write. Click View to open a closable preview window.")
    for name in order:
        row = st.columns([0.9, 3.1], gap="small")
        if row[0].button("View", key=f"view_{name}"):
            show_preview(name, files[name])
        row[1].write(f"`{name}`")

elif page == "Star":
    st.subheader("Star")
    star = system["star"]
    st.caption("Fit flags control the first two columns of the star.dat row in bodies.lst. Prior flags add entries to priors.in.")
    header = st.columns([1.5, 2, 2, 1.5, 1.5])
    for col, text in zip(header, ["Parameter", "Value", "Sigma", "Fit", "Prior"]):
        col.markdown(f"**{text}**")
    for key, title, unit in [("mass", "Mass", "Msun"), ("radius", "Radius", "Rsun")]:
        row = st.columns([1.5, 2, 2, 1.5, 1.5])
        row[0].write(f"{title} [{unit}]")
        with row[1]:
            value_field = float_input("", star[key], f"star_{key}_value")
        with row[2]:
            sigma_field = float_input("", star[f"{key}_err"], f"star_{key}_sigma", min_value=0.0)
        with row[3]:
            fit_field = st.checkbox("Fit", value=bool(star.get(f"fit_{key}", False)), key=f"star_{key}_fit", label_visibility="collapsed")
        with row[4]:
            prior_field = st.checkbox("Prior", value=bool(star.get(f"{key}_prior", False)), key=f"star_{key}_prior", label_visibility="collapsed")
        star[key] = value_field
        star[f"{key}_err"] = sigma_field
        star[f"fit_{key}"] = fit_field
        star[f"{key}_prior"] = prior_field

elif page == "Planets":
    st.subheader("Planets")
    st.caption("Names are lower-case by design and default to b, c, d, ... . Fit and prior controls follow the bodies.lst parameter order.")
    ensure_planets(system)
    for index, planet in enumerate(system["planets"]):
        default_label = get_planet_label(index)
        with st.expander(f"Planet {index + 1}", expanded=index == 0):
            label = st.text_input("Planet name", value=planet.get("label", default_label), key=f"label_{index}")
            planet["label"] = label.lower() or default_label
            planet["transiting"] = st.checkbox("Transiting", value=planet.get("transiting", True), key=f"transit_{index}")
            parameter_tab, transit_tab = st.tabs(["Parameters", "Transit times"])
            with parameter_tab:
                header = st.columns([1.5, 1.3, 1.3, 1.3, 0.7, 0.7, 1.3, 1.3])
                for col, text in zip(header, ["Parameter", "Value", "Min", "Max", "Fit", "Prior", "Prior value", "Prior sigma"]):
                    col.markdown(f"**{text}**")
                parameter_titles = {
                    "mass": f"Mass [m{index + 2} (Mjup)]",
                    "radius": f"Radius [R{index + 2} (Rjup)]",
                    "period": f"Period [P (days)]",
                    "ecc": f"Eccentricity [e{index + 2}]",
                    "omega": f"Argument of pericentre [w{index + 2} (deg)]",
                    "M0": f"Mean anomaly [mA{index + 2} (deg)]",
                    "inc": f"Inclination [i{index + 2} (deg)]",
                    "lnode": f"Longitude of ascending node [lN{index + 2} (deg)]",
                }
                for key, prior_label, title in FIT_PARAMETERS:
                    row = st.columns([1.5, 1.3, 1.3, 1.3, 0.7, 0.7, 1.3, 1.3])
                    row[0].write(parameter_titles[key])
                    with row[1]:
                        value_field = float_input("", planet[key], f"{key}_value_{index}")
                    with row[2]:
                        min_field = float_input("", planet[f"{key}_min"], f"{key}_min_{index}")
                    with row[3]:
                        max_field = float_input("", planet[f"{key}_max"], f"{key}_max_{index}")
                    with row[4]:
                        fit_field = st.checkbox("Fit", value=bool(planet.get(f"fit_{key}", True)), key=f"{key}_fit_{index}", label_visibility="collapsed")
                    with row[5]:
                        prior_field = st.checkbox("Prior", value=bool(planet.get(f"prior_{key}", False)), key=f"{key}_prior_{index}", label_visibility="collapsed")
                    with row[6]:
                        prior_value_field = float_input("", planet.get(f"prior_{key}_value", planet[key]), f"{key}_prior_value_{index}")
                    with row[7]:
                        sigma_field = float_input("", planet.get(f"prior_{key}_sigma", 0.0), f"{key}_sigma_{index}", min_value=0.0)
                    planet[key] = value_field
                    planet[f"{key}_min"] = min_field
                    planet[f"{key}_max"] = max_field
                    planet[f"fit_{key}"] = fit_field
                    planet[f"prior_{key}"] = prior_field
                    planet[f"prior_{key}_value"] = prior_value_field
                    planet[f"prior_{key}_sigma"] = sigma_field
                st.caption("Semi-major axis and time of pericentre passage remain derived .dat rows and are not fit columns in bodies.lst.")
            with transit_tab:
                edit_transit_times(planet, index + 2)

elif page == "RV Observations":
    st.subheader("RV Observations")
    observations = system.get("rv_obs", [])
    count = st.number_input("Number of RV measurements", 0, 1000, len(observations), key="rv_count")
    rv_headers = ["Time_RV(days)", "RV(m/s)", "err_RV(m/s)", "RV_set"]
    header_cols = st.columns(4)
    for column, header in zip(header_cols, rv_headers):
        column.markdown(f"**{header}**")
    updated = []
    for index in range(count):
        old = observations[index] if index < len(observations) else {}
        cols = st.columns(4)
        with cols[0]:
            jd_field = float_input("", old.get("jd", 2455000.0), f"rv_jd_{index}")
        with cols[1]:
            rv_field = float_input("", old.get("rv", 0.0), f"rv_{index}")
        with cols[2]:
            erv_field = float_input("", old.get("erv", 1.0), f"erv_{index}")
        rvset_field = cols[3].number_input("RV set", value=int(old.get("rvset", 1)), key=f"rvset_{index}", label_visibility="collapsed")
        updated.append({"jd": jd_field, "rv": rv_field, "erv": erv_field, "rvset": rvset_field})
    system["rv_obs"] = updated

elif page == "Integration Settings":
    st.subheader("Integration Settings (arg.in)")
    args = system.setdefault("arg_in", {})
    ensure_arg_defaults(args)
    header = st.columns([1.5, 1.5, 5])
    for col, text in zip(header, ["Parameter", "Value", "Description"]):
        col.markdown(f"**{text}**")
    for key, ptype, default, comment, choices in ARG_PARAMETERS:
        skey = f"arg_{key}"
        row = st.columns([1.5, 1.5, 5])
        row[0].write(key)
        current = args.get(key, default)
        with row[1]:
            if ptype == "bool":
                value = st.checkbox("", value=bool(current), key=skey, label_visibility="collapsed")
            elif ptype == "choice":
                index = choices.index(current) if current in choices else 0
                value = st.selectbox("", choices, index, key=skey, label_visibility="collapsed")
            elif ptype == "int":
                value = int(st.number_input("", value=int(current), key=skey, label_visibility="collapsed"))
            else:
                value = float_input("", float(current), skey)
        row[2].caption(comment)
        args[key] = value
    nb = system["n_planets"] + 1
    row = st.columns([1.5, 1.5, 5])
    row[0].write("NB")
    row[1].number_input("", value=nb, key="arg_nb", disabled=True, label_visibility="collapsed")
    row[2].caption(NB_COMMENT)
    args["NB"] = nb

elif page == "Analysis Settings":
    st.subheader("Analysis Settings (configuration.yml)")
    cfg = system.setdefault("config", {})
    ensure_config_defaults(cfg)
    tab_run, tab_ana, tab_oc, tab_rv = st.tabs(["run", "analysis", "OC", "RV"])
    with tab_run:
        render_params(cfg, "run", RUN_PARAMS)
        st.subheader("PyDE")
        render_params(cfg, "pyde", PYDE_PARAMS)
        st.subheader("emcee")
        render_params(cfg, "emcee", EMCEE_PARAMS)
    with tab_ana:
        render_params(cfg, "analysis", ANALYSIS_PARAMS)
    with tab_oc:
        render_params(cfg, "oc", OC_PARAMS)
    with tab_rv:
        render_params(cfg, "rv", RV_PARAMS)
