#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 18:03:09 2025

@author: harshada
"""
# edge_interactions.py
import pandas as pd
import os

def process_edge_directions(input_file, output_dir):
    df = pd.read_csv(input_file, sep="\t")

    def process_direction(row):
        source, target, annotation, direction = row["source"], row["target"], row["Annotation"], row["Direction"]
        results = []
        annotation_lower = annotation.strip().lower()

        if annotation_lower in ["complex", "complex; input", "reaction", "complex; input; reaction"]:
            results.append({"source": source, "target": target, "interaction": "->", "Annotation": annotation + " (Complex interaction)", "Main Direction": "->"})
            return results
        if "PPrel: binding/association" in annotation and direction == "-":
            results.append({"source": source, "target": target, "interaction": "->", "Annotation": annotation + " (Binding/Association)", "Main Direction": "->"})
            return results
        if annotation_lower in ["input", "ecrel: compound", "pprel: compound", "complex; reaction", "pprel"] and direction == "-":
            results.append({"source": source, "target": target, "interaction": "->", "Annotation": annotation + " (Input/Compound interaction)", "Main Direction": "->"})
            return results

        if direction == "<-":
            results.append({"source": target, "target": source, "interaction": "->", "Annotation": annotation + " (Swapped)", "Main Direction": "->"})
        elif direction == "<-|":
            results.extend([
                {"source": target, "target": source, "interaction": "->", "Annotation": annotation + " (Target activates Source)", "Main Direction": "->"},
                {"source": source, "target": target, "interaction": "-|", "Annotation": annotation + " (Source inhibits Target)", "Main Direction": "-|"},
            ])
        elif direction == "<->":
            results.extend([
                {"source": source, "target": target, "interaction": "->", "Annotation": annotation + " (Bidirectional Activation)", "Main Direction": "->"},
                {"source": target, "target": source, "interaction": "->", "Annotation": annotation + " (Bidirectional Activation)", "Main Direction": "->"},
            ])
        elif direction == "|->":
            results.extend([
                {"source": source, "target": target, "interaction": "->", "Annotation": annotation + " (Source activates Target)", "Main Direction": "->"},
                {"source": target, "target": source, "interaction": "-|", "Annotation": annotation + " (Target inhibits Source)", "Main Direction": "-|"},
            ])
        elif direction == "|-":
            results.append({"source": target, "target": source, "interaction": "-|", "Annotation": annotation + " (Target inhibits Source)", "Main Direction": "-|"})
        elif direction == "-|":
            results.append({"source": source, "target": target, "interaction": "-|", "Annotation": annotation, "Main Direction": "-|"})
        else:
            results.append({"source": source, "target": target, "interaction": direction, "Annotation": annotation, "Main Direction": direction})

        return results

    processed_rows = [entry for _, row in df.iterrows() for entry in process_direction(row)]
    df_final = pd.DataFrame(processed_rows)
    df_final = df_final[df_final["interaction"] != "-"]

    os.makedirs(output_dir, exist_ok=True)

    # Add regulation
    df_final["Regulation"] = df_final["interaction"].map({"->": "Activation", "-|": "Inhibition"})
    df_final.drop(columns="Main Direction", inplace=True)

    final_output_file = os.path.join(output_dir, "finalized_interactions_with_regulation.tsv")
    df_final.to_csv(final_output_file, sep="\t", index=False)

    return final_output_file
