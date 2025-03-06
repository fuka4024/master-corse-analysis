import os
from pymol import cmd
import csv
import math

base_dir = "C:/Users/monom/Downloads/BmToll_dim_spz3_dim/folds_2024_12_16_12_00/"
output_dir = "C:/Users/monom/Documents/Toll_pre/"

# フォルダ内のすべてのCIFファイルを処理する
for folder_name in os.listdir(base_dir):
    folder_path = os.path.join(base_dir, folder_name)
    if os.path.isdir(folder_path):  # フォルダかどうか確認
        for file_name in os.listdir(folder_path):
            if file_name.endswith(".cif"):  # CIFファイルのみ処理
                cif_path = os.path.join(folder_path, file_name)
                base_name = os.path.splitext(file_name)[0]  # ファイル名（拡張子なし）
                
                # PyMOLの初期化
                cmd.reinitialize()
                cmd.load(cif_path, base_name)

                # タンパク質のロードと選択
                cmd.select("Spz", "chain B or chain C")
                cmd.select("Toll", "chain A")

                # 最短距離
                min_distance = cmd.distance("closest_atoms", "Spz", "Toll", 4.0)

                # 平均距離計算
                def calculate_average_distance(selection1, selection2, cutoff=4.0):
                    distances = []
                    model1 = cmd.get_model(selection1)
                    model2 = cmd.get_model(selection2)
                    for atom1 in model1.atom:
                        for atom2 in model2.atom:
                            dist = math.sqrt(
                                (atom1.coord[0] - atom2.coord[0])**2 +
                                (atom1.coord[1] - atom2.coord[1])**2 +
                                (atom1.coord[2] - atom2.coord[2])**2
                            )
                            if dist <= cutoff:
                                distances.append(dist)
                    return sum(distances) / len(distances) if distances else None

                average_distance = calculate_average_distance("Spz", "Toll", cutoff=4.0)

                # 重心間距離
                def calculate_com(selection):
                    model = cmd.get_model(selection)
                    x, y, z, mass_sum = 0, 0, 0, 0
                    for atom in model.atom:
                        mass = atom.get_mass()
                        x += atom.coord[0] * mass
                        y += atom.coord[1] * mass
                        z += atom.coord[2] * mass
                        mass_sum += mass
                    return (x / mass_sum, y / mass_sum, z / mass_sum)

                com_Spz = calculate_com("Spz")
                com_Toll = calculate_com("Toll")
                com_distance = math.sqrt(
                    (com_Spz[0] - com_Toll[0])**2 +
                    (com_Spz[1] - com_Toll[1])**2 +
                    (com_Spz[2] - com_Toll[2])**2
                )

                # Interface面積の計算
                cmd.select("interface_Spz", "Spz within 5 of Toll and (Spz and not Toll)")
                cmd.select("interface_Toll", "Toll within 5 of Spz and (Toll and not Spz)")

                # Spzの面積を計算
                area_Spz = cmd.get_area("interface_Spz and Spz")
                # Tollの面積を計算
                area_Toll = cmd.get_area("interface_Toll and Toll")

                # 合計Interface面積
                interface_area = area_Spz + area_Toll

                # CSV出力
                output_path_csv = os.path.join(output_dir, f"{folder_name}_{base_name}_binding_sites.csv")
                with open(output_path_csv, "w", encoding="utf-8", newline="") as f:
                    writer = csv.writer(f)
                    writer.writerow(["Metric", "Value"])
                    writer.writerow(["Shortest Distance (Å)", min_distance])
                    writer.writerow(["Average Distance (Å)", average_distance])
                    writer.writerow(["Center of Mass Distance (Å)", com_distance])
                    writer.writerow(["Interface Area (Å²)", interface_area])

                # PyMOL画像保存
                def highlight_binding_sites(selection1, selection2, cutoff=4.0):
                    cmd.select("binding_sites", f"({selection1}) within {cutoff} of ({selection2})")  # 範囲内を選択

                # 結合箇所を表示 (Spz のみ強調表示)
                highlight_binding_sites("Toll", "Spz", cutoff=4.0)

                cmd.color("blue", "Spz")
                cmd.color("green", "Toll")
                cmd.show("cartoon", "Spz or Toll")
                cmd.show("spheres", "binding_sites")  # スフィアで表示
                cmd.color("orange", "binding_sites")
                cmd.hide("labels")
                cmd.hide("dashes", "closest_atoms")
                cmd.bg_color("white")
                cmd.ray()

                # 出力画像のパス
                output_path_png = os.path.join(output_dir, f"{folder_name}_{base_name}_binding_sites.png")
                cmd.png(output_path_png)

                print(f"Processed {file_name} in folder {folder_name} and saved images.")
