import json
import pytest
from qstars.cli import main
from qstars.config import resolve
from qstars.evolution import run

def test_reference_cli(tmp_path):
    assert main(['reference','benchmarks/ball2011_polytrope.yaml','--output',str(tmp_path)])==0
    d=next(tmp_path.iterdir())
    assert (d/'sequence.csv').is_file()
    assert json.loads((d/'termination.json').read_text())['scope'].startswith('Dimensionless')

def test_invalid_physics_not_silently_accepted():
    with pytest.raises(ValueError):resolve({'boundary':{'inner':'full_saturated_model'}})
    with pytest.raises(ValueError):resolve({'evolution':{'max_wall_seconds':{'kind':'pulse'}}})

def test_wall_budget_preserves_restart(tmp_path):
    p,r=run(resolve({'evolution':{'max_wall_seconds':1e-12}}),tmp_path)
    assert r['status']=='control_stop'
    assert (p/'checkpoint.npz').exists()
