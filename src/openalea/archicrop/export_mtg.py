""" This module provides functions to save MTG and PlantGL scene into MTG and OBJ files.

```python
    fn_mtg = Path(output_fn).with_suffix('.mtg')
    fn_obj = Path(output_fn).with_suffix('.obj')
    save_mtg(g, scene, fn_mtg, fn_obj)
```

"""
from pathlib import Path
from openalea.mtg.io import write_mtg
from openalea.plantgl.all import Scene


def fix_obj(fn):
    """Fix the OBJ file by removing 'vn ' and rewriting 'f ' without 'vn '."""
    with open(fn, "r") as f:
        lines = f.readlines()

    out = []

    for line in lines:
        if line.startswith("vn "):
            continue
        if line.startswith("f "):
            tokens = line.split()
            new_face = ["f"]

            for t in tokens[1:]:
                parts = t.split("/")
                if len(parts) == 3:
                    # v/vt/vn -> v/vt
                    if parts[1]:
                        new_face.append("/".join(parts[:2]))
                    else:
                        new_face.append(parts[0])
                elif len(parts) == 2:
                    new_face.append(t)
                else:
                    new_face.append(parts[0])

            out.append(" ".join(new_face) + "\n")
        else:
            sline = line.strip()
            if sline:
                out.append(line)

    with open(fn, "w") as f:
        f.writelines(out)


def save_mtg(g, scene, mtg_fn, obj_fn):
    """Saves MTG and PlantGL scene into MTG and OBJ files."""
    g.properties()['Id']={}
    ids = g.properties()['Id']
    for vid in g.vertices(): 
        ids[vid]=vid

    def bools(g):
        for p in ['is_green', 'grow', 'dead', 'senescence']:
            _prop = g.properties()[p]
            for vid in _prop:
                _prop[vid] = int(bool(_prop[vid]))
    
    bools(g)
    props = [
        ('Id','INT'),
     ('rank', 'INT'),
     ('length','REAL'),
     ('visible_length','REAL'),
     ('is_green','INT'),
     ('stem_diameter','REAL'),
     ('azimuth','REAL'),
     ('grow','INT'),
     ('age','REAL'),
     ('tiller_angle','REAL'),
     ('leaf_area','REAL'),
     ('visible_leaf_area','REAL'),
     ('senescent_area','REAL'),
     ('senescent_length','REAL'),
     ('shape_max_width','REAL'),
     ('dead','INT'),
     ('inclination','REAL'),
     ('senescence', 'INT')
    ] 

    mtg_lines=write_mtg(g, props)
    fn = mtg_fn
    if fn.exists():
        fn.unlink(missing_ok=True)

    fn.write_text(mtg_lines)
    print(f'Write MTG in {str(fn)}')

    senescent_scene = Scene()
    # use scene.todict() rather that iterating on the scene
    for vid, shs in scene.todict().items():
        if len(shs) == 2:
            # check is senecescent
            sh, sh_sen = shs
            sh_sen.setName(f'sen_vid_{vid}')
            print(f'vid {vid} has 2 shapes, senescent and non-senescent')
        else:
            sh = shs[0]           
            
        sh.setName(f'vid_{vid}')

    fn = obj_fn
    if fn.exists():
        fn.unlink(missing_ok=True)
    scene.save(str(obj_fn))

    fix_obj(obj_fn)
    
