import glob
import os
from qgis.core import QgsVectorLayer, QgsProject

folder = '/Volumes/Sandisk4TB/PhD_MS/TARSAN/Radar/GSSI/2022_2023/raw/KML_THW/'   # <-- change this

for csv_path in glob.glob(os.path.join(folder, "*.csv")):
    uri = f"file:///{csv_path}?delimiter=,&xField=lon&yField=lat&crs=EPSG:4326"
    layer_name = os.path.basename(csv_path)
    layer = QgsVectorLayer(uri, layer_name, "delimitedtext")

    if layer.isValid():
        QgsProject.instance().addMapLayer(layer)
        print("Loaded:", layer_name)
    else:
        print("FAILED:", layer_name)
