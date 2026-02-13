lmp_new -in conv_cfg_script
python free.py --use-radii --radii "Cu=1.28,Al=1.43,Ti=1.47,Zr=1.60,Ag=1.44" --export-per-atom final.cfg
python rdf.py
python bond.py final.cfg
python voronoi.py
