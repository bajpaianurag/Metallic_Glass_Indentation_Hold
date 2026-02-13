python d2min.py
python load_depth_headness_modulus.py > indent_summary
cp ../final.cfg final.cfg
python bond.py final.cfg
python voronoi.py
