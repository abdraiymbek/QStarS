"""Serializable controls; no eval or executable configuration."""
import numpy as np

def value(control, t_years):
    if isinstance(control,(int,float)):
        return float(control)
    kind = control["kind"]
    if kind == "table":
        return float(np.interp(t_years,control["time_years"],control["values"]))
    if kind == "pulse":
        return float(control["base"] + (control["amplitude"] if
            control["start_years"] <= t_years < control["end_years"] else 0.))
    if kind == "sinusoid":
        return float(control["base"]+control["amplitude"]*np.sin(
            2*np.pi*t_years/control["period_years"]+control.get("phase",0)))
    raise ValueError(f"Unknown control {kind}")

def next_event(control,t):
    if not isinstance(control,dict): return float("inf")
    points = control.get("time_years",[control.get("start_years",float("inf")),
                                     control.get("end_years",float("inf"))])
    return min((x for x in points if x > t+1e-13), default=float("inf"))
