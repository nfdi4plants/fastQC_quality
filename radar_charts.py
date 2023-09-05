#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 08:30:09 2023

@author: hannah
"""

import zipfile
import os
from csv import reader
import pandas as pd
import plotly.graph_objects as go
# import plotly.offline as pyo
import plotly.io as io
import numpy as np
import matplotlib.pyplot as plt
import glob
import shutil
import shapely.geometry


def getPASScolors(df):
    colors = []
    for elem in df["summaries"]:
        if elem == "PASS":
            colors.append("limegreen")
        elif elem == "WARN":
            colors.append("orange")
        else:
            colors.append("red")
    colors.append(colors[0])
    return colors

def readFile(direct, file):
    temp = direct + "tmp"
    isExists = os.path.exists(temp)
    if not isExists:
        os.makedirs(temp)
    filename = file[8:-4]
    # print(direct)
    infoList = []
    baseStatistics = None
    perBaseSequenceQuality = None
    perTileSequenceQuality = None
    perSequenceQualityScore = None
    perBaseSequenceContent = None
    perSequenceGCContent = None
    perBaseNContent = None
    sequenceLengthDistribution = None
    with zipfile.ZipFile(file) as zip_ref:
        zip_ref.extractall(temp)

    # print(temp, " ", filename, "\n")
    with open(temp+filename+"/fastqc_data.txt", "r") as f:
          file_reader = reader(f, delimiter="\t")
          switch = ""
          for row in file_reader:
              if row[0] == ">>Basic Statistics":
                  switch = "baseStatistics"
                  infoList = []
              if row[0] == ">>Per base sequence quality":
                  switch = "perBaseSequenceQuality"
                  infoList = []
              if row[0] == ">>Per tile sequence quality":
                  switch = "perTileSequenceQuality"
                  infoList = []
              if row[0] == ">>Per sequence quality scores":
                  switch = "perSequenceQualityScore"
                  infoList = []
              if row[0] == ">>Per base sequence content":
                  switch = "perBaseSequenceContent"
                  infoList = []
              if row[0] == ">>Per sequence GC content":
                  switch = "perSequenceGCContent"
                  infoList = []
              if row[0] == ">>Per base N content":
                  switch = "perBaseNContent"
                  infoList = []
              if row[0] == ">>Sequence Length Distribution":
                  switch = "sequenceLengthDistribution"
                  infoList = [] 
                  
              if row[0] == ">>END_MODULE":
                  if switch == "baseStatistics":
                      colnames = infoList[0]
                      infoList.pop(0)
                      baseStatistics = pd.DataFrame(infoList, columns=colnames)
                  if switch == "perBaseSequenceQuality":
                      colnames = infoList[0]
                      infoList.pop(0)
                      perBaseSequenceQuality = pd.DataFrame(infoList, columns=colnames)
                  if switch == "perTileSequenceQuality":
                      colnames = infoList[0]
                      infoList.pop(0)
                      perTileSequenceQuality = pd.DataFrame(infoList, columns=colnames)
                  if switch == "perSequenceQualityScore":
                      colnames = infoList[0]
                      infoList.pop(0)
                      perSequenceQualityScore = pd.DataFrame(infoList, columns=colnames)
                  if switch == "perBaseSequenceContent":
                      colnames = infoList[0]
                      infoList.pop(0)
                      perBaseSequenceContent = pd.DataFrame(infoList, columns=colnames)
                  if switch == "perSequenceGCContent":
                      colnames = infoList[0]
                      infoList.pop(0)
                      perSequenceGCContent = pd.DataFrame(infoList, columns=colnames)
                  if switch == "perBaseNContent":
                      colnames = infoList[0]
                      infoList.pop(0)
                      perBaseNContent = pd.DataFrame(infoList, columns=colnames)
                  if switch == "sequenceLengthDistribution":
                      colnames = infoList[0]
                      infoList.pop(0)
                      sequenceLengthDistribution = pd.DataFrame(infoList, columns=colnames) 
                  switch = "endModule"
    
              if switch == "baseStatistics":
                  if row[0] == ">>Basic Statistics":
                      pass
                  else:
                      infoList.append(row)
              if switch == "perBaseSequenceQuality":
                  if row[0] == ">>Per base sequence quality":
                      pass
                  elif row[0] == "#Base":
                      infoList.append(row)
                  elif "-" in row[0]:
                      bases = row[0].split("-")
                      bases=range(int(bases[0]), int(bases[1])+1)
                      i=0
                      for elem in bases:
                          row = [float(x) for x in row[1:]]
                          row.insert(0, float(bases[i]))
                          infoList.append(row)
                          i+=1
                  else:
                      row = [float(x) for x in row]
                      infoList.append(row)
              if switch == "perTileSequenceQuality":
                  if row[0] == ">>Per tile sequence quality":
                      pass
                  elif row[0] == "#Tile":
                      infoList.append(row)
                  elif "-" in row[1]:
                      bases = row[1].split("-")
                      bases=range(int(bases[0]), int(bases[1])+1)
                      tile=float(row[0])
                      mean=float(row[2])
                      i=0
                      for elem in bases:
                          row = [tile, float(bases[i]) , mean]
                          # print(row)
                          infoList.append(row)
                          i+=1
                  else:
                      row = [float(x) for x in row]
                      infoList.append(row)
              if switch == "perSequenceQualityScore":
                  if row[0] == ">>Per sequence quality scores":
                      pass
                  elif row[0] == "#Quality":
                      infoList.append(row)
                  else:
                      row = [float(x) for x in row]
                      infoList.append(row)
              if switch == "perBaseSequenceContent":
                  if row[0] == ">>Per base sequence content":
                      pass
                  elif row[0] == "#Base":
                      infoList.append(row)
                  elif "-" in row[0]:
                      bases = row[0].split("-")
                      bases=range(int(bases[0]), int(bases[1])+1)
                      i=0
                      for elem in bases:
                          row = [float(x) for x in row[1:]]
                          row.insert(0, float(bases[i]))
                          infoList.append(row)
                          i+=1
                  else:
                      row = [float(x) for x in row]
                      infoList.append(row)
              if switch == "perSequenceGCContent":
                  if row[0] == ">>Per sequence GC content":
                      pass
                  elif row[0] == "#GC Content":
                      infoList.append(row)
                  else:
                      row = [float(x) for x in row]
                      infoList.append(row)
              if switch == "perBaseNContent":
                  if row[0] == ">>Per base N content":
                      pass
                  elif row[0] == "#Base":
                      infoList.append(row)
                  elif "-" in row[0]:
                      bases = row[0].split("-")
                      bases=range(int(bases[0]), int(bases[1])+1)
                      i=0
                      for elem in bases:
                          row = [float(x) for x in row[1:]]
                          row.insert(0, float(bases[i]))
                          infoList.append(row)
                          i+=1
                  else:
                      row = [float(x) for x in row]
                      infoList.append(row)
              if switch == "sequenceLengthDistribution":
                  if row[0] == ">>Sequence Length Distribution":
                      pass
                  elif row[0] == "#Length":
                      infoList.append(row)
                  elif "-" in row[0]:
                      bases = row[0].split("-")
                      bases=range(int(bases[0]), int(bases[1])+1)
                      i=0
                      for elem in bases:
                          if i == 0:
                              row = [float(x) for x in row[1:]]
                              row.insert(0, int(bases[0]))
                          else:
                              row = [elem, float(0.0)]
                          infoList.append(row)
                          i+=1
                  else:
                      row = [float(x) for x in row]
                      infoList.append(row) 
                      
    summary = pd.read_csv(temp+filename+"/summary.txt", sep="\t")
    
    shutil.rmtree(temp)
    return baseStatistics, perBaseSequenceQuality, perTileSequenceQuality, perSequenceQualityScore, perBaseSequenceContent, perSequenceGCContent, perBaseNContent, sequenceLengthDistribution, summary


def scorePerBaseSequenceQuality(df, seq_length):
    scores = []
    for elem in df["Median"]:
        if elem < 20: # fail
            scores.append(0)
        elif elem < 25: # warn
            scores.append(0.5)
        else: # pass
            scores.append(1)
    i=0
    for elem in df["Lower Quartile"]:
        if elem < 5: # fail
            scores[i] = 0
        elif elem < 10 and scores[i] < 0: # warn
            scores[i] = 0.5
        i += 1
    qs = (sum(scores)/seq_length)*100
    return(qs)

def scorePerTileSequenceQuality(df, seq_length):
    means = []
    score = 0
    for base in range(seq_length):
        baseMean = df.loc[df["Base"] == base+1, "Mean"].sum()/(len(df)/seq_length)
        means.append(baseMean)
    i=0
    for elem in df["Mean"]:
        if elem < means[i]-5: # fail
            pass
        elif elem < means[i]-2: # warn
            score += 0.5
        else: # pass
            score += 1
        if i == 7:
            i = 0
        else:
            i += 1
    qs = (score/len(df))*100
    return qs


def scorePerSequenceQualityScore(df, tot_seqs):
    # average = max(df["Count"])
    # quality = df.loc[df["Count"] == average, "#Quality"].iloc[0]
    sums = 0
    for elem in df["#Quality"]:
        if elem < 20: # fail
            pass
        elif elem < 27: # warn
            sums = sums + (df.loc[df["#Quality"] == elem, "Count"].iloc[0] * 0.5)
        else: # pass
            sums = sums + df.loc[df["#Quality"] == elem, "Count"].iloc[0]
    qs = (sums/tot_seqs)*100
    return qs
    

def scorePerBaseSequenceContent(df, seq_length):
    score = 0
    for elem in df["#Base"]:
        at_score = abs(df.loc[df["#Base"] == elem, "A"].iloc[0] - df.loc[df["#Base"] == elem, "T"].iloc[0])
        gc_score = abs(df.loc[df["#Base"] == elem, "G"].iloc[0] - df.loc[df["#Base"] == elem, "C"].iloc[0])
        if at_score > 20 or gc_score > 20: # fail
            pass
        elif at_score > 10 or gc_score > 10: # warn
            score += 0.5
        else: # pass
            score += 1
    qs = (score/seq_length)*100
    return qs


def scorePerBaseNContent(df, seq_length):
    sums = 0
    for elem in df["N-Count"]:
        if elem > 20: # fail
            pass
        elif elem > 5: # warn
            sums += 0.5
        else: # pass
            sums += 1
    qs = (sums/seq_length)*100
    return qs


def scoreSequenceLengthDistribution(df, tot_seqs):
    average = max(df["Count"])
    length_sum = 0
    # quality = df.loc[df["Count"] == average, "#Quality"].iloc[0]

    for elem in df["#Length"]:
        length_sum += (elem*df.loc[df["#Length"] == elem, "Count"].iloc[0])
    real_mean = length_sum/tot_seqs
    qs = (average/tot_seqs)*100
    return qs


def makeRadarChart(df):
    categories = list(df["field"])
    categories.append(df["field"].loc[0])
    qualities = list(df["qs"])
    qualities.append(df["qs"].loc[0])
    label_loc = np.linspace(start=0, stop=2 * np.pi, num=len(qualities))
    
    plt.figure(figsize=(8,8))
    plt.subplot(polar=True)
    plt.plot(label_loc, qualities, label="test")
    lines, labels = plt.thetagrids(np.degrees(label_loc), labels=categories)
    plt.legend()
    plt.show()

def makeOtherRadarChart(df, file, colors):
    filename = file[9:-11]
    minimum = round(min(df["qs"])-(min(df["qs"])%10))
    # print(minimum)
    categories = list(df["field"])
    categories.append(df["field"].loc[0])
    qualities = list(df["qs"])
    qualities.append(df["qs"].loc[0])
    # print(colors)
    
    fig = go.Figure(
        data=[
            go.Scatterpolar(r=qualities, theta=categories, fill="toself", name="test", 
                            marker=dict(color=colors, size=10))
            ],
        layout=go.Layout(
            title=go.layout.Title(text=filename),
            polar={"radialaxis": {"visible": True,
                                  "range": (minimum,100)}},
            showlegend=False
            )
        )
    static_image_bytes = fig.to_image(format="png", scale = 4, engine="kaleido")
    # fig.show()
    # pyo.plot(fig)
    #save_loc = file[:-4].split("r_fastqc/")[0] + "r_quality/" + file[:-4].split("r_fastqc/")[1]
    save_loc = "./outputs/" + file[9:-4]
    print(save_loc)
    # fig.write_image(save_loc+"_radarplot.png")
    with open(save_loc+"_radarplot.png", "wb") as f:
        f.write(static_image_bytes)

def calc_score(df, file):
    df2 = df.drop(df[df["field"] == "Per Tile Sequence Quality"].index)
    # print(df2)
    df2["trace"] = [1,1,1,1,1]
    df2["theta_n"] = pd.factorize(df2["field"])[0]
    df2["theta_radian"] = (df2["theta_n"] / (df2["theta_n"].max() + 1)) * 2 * np.pi
    df2["x"] = np.cos(df2["theta_radian"]) * df2["qs"]
    df2["y"] = np.sin(df2["theta_radian"]) * df2["qs"]
    # print(df)
    df_a = df2.groupby("trace").apply(
        lambda d: shapely.geometry.MultiPoint(list(zip(d["x"], d["y"]))).convex_hull.area
        )
    filename = file[:-11]
    # save_loc = file[:-4].split("r_fastqc/")[0] + "r_quality/" + file[:-4].split("r_fastqc/")[1]
    save_loc = "./outputs/"+filename[9:]
    score = str(round(df_a.values[0],2))
    print(save_loc, "_area.txt")
    with open(save_loc+"_area.txt", "w") as f:
        f.write(score)
    


directory = "./inputs/"
# a=os.listdir(directory)
# print(a)
# print(directory)
for file_path in glob.glob(os.path.join(directory, "*zip")):
    print(file_path)
    baseStatistics, perBaseSequenceQuality, perTileSequenceQuality, perSequenceQualityScore, perBaseSequenceContent, perSequenceGCContent, perBaseNContent, sequenceLengthDistribution, summary = readFile(directory, file_path)
    t_total_seqs = baseStatistics.loc[baseStatistics["#Measure"] == "Total Sequences", "Value"]
    total_seqs = int(t_total_seqs.iloc[0])
    if "-" in baseStatistics.loc[baseStatistics["#Measure"] == "Sequence length", "Value"].iloc[0]:
        print(f"WARNING: All sequences are not of same length! The troublesome file is: {file_path}")
        continue
    t_sequence_length = baseStatistics.loc[baseStatistics["#Measure"] == "Sequence length", "Value"]
    sequence_length = int(t_sequence_length.iloc[0])
    qs_2 = scorePerBaseSequenceQuality(df=perBaseSequenceQuality, seq_length=sequence_length)
    summary_2 = summary.loc[summary["Basic Statistics"] == "Per base sequence quality", "PASS"].iloc[0]
    if perTileSequenceQuality is not None:
        qs_3 = scorePerTileSequenceQuality(df=perTileSequenceQuality, seq_length=sequence_length)
        summary_3 = summary.loc[summary["Basic Statistics"] == "Per tile sequence quality", "PASS"].iloc[0]
    qs_4 = scorePerSequenceQualityScore(df=perSequenceQualityScore, tot_seqs=total_seqs)
    summary_4 = summary.loc[summary["Basic Statistics"] == "Per sequence quality scores", "PASS"].iloc[0]
    qs_5 = scorePerBaseSequenceContent(df=perBaseSequenceContent, seq_length=sequence_length)
    summary_5 = summary.loc[summary["Basic Statistics"] == "Per base sequence content", "PASS"].iloc[0]
    # scorePerSequenceGCContent(df=perSequenceGCContent)
    qs_7 = scorePerBaseNContent(df=perBaseNContent, seq_length=sequence_length)
    summary_7 = summary.loc[summary["Basic Statistics"] == "Per base N content", "PASS"].iloc[0]
    qs_8 = scoreSequenceLengthDistribution(df=sequenceLengthDistribution, tot_seqs=total_seqs)
    summary_8 = summary.loc[summary["Basic Statistics"] == "Sequence Length Distribution", "PASS"].iloc[0]
    if perTileSequenceQuality is not None:
        qs_df = pd.DataFrame(dict(
            qs=[qs_2, qs_4, qs_5, qs_7, qs_8], 
            # qs=[qs_2, qs_3, qs_4, qs_5, qs_7, qs_8], 
            # field=["Per Base Sequence Quality", "Per Tile Sequence Quality", "Per Sequence Quality Score", "Per Base Sequence Content", "Per Base N content", "Sequence Length Distribution"]))
            field=["Per Base Sequence Quality", "Per Sequence Quality Score", "Per Base Sequence Content", "Per Base N content", "Sequence Length Distribution"]))
        summary_df = pd.DataFrame(dict(
            summaries=[summary_2, summary_4, summary_5, summary_7, summary_8],
            # summaries=[summary_2, summary_3, summary_4, summary_5, summary_7, summary_8],
            # field=["Per Base Sequence Quality", "Per Tile Sequence Quality", "Per Sequence Quality Score", "Per Base Sequence Content", "Per Base N content", "Sequence Length Distribution"]))
            field=["Per Base Sequence Quality", "Per Sequence Quality Score", "Per Base Sequence Content", "Per Base N content", "Sequence Length Distribution"]))
    else:
        qs_df = pd.DataFrame(dict(
            qs=[qs_2, qs_4, qs_5, qs_7, qs_8], 
            field=["Per Base Sequence Quality", "Per Sequence Quality Score", "Per Base Sequence Content", "Per Base N content", "Sequence Length Distribution"]))
        summary_df = pd.DataFrame(dict(
            summaries=[summary_2, summary_4, summary_5, summary_7, summary_8],
            field=["Per Base Sequence Quality", "Per Sequence Quality Score", "Per Base Sequence Content", "Per Base N content", "Sequence Length Distribution"]))
    # # qs_df["qs"] = [0, 0, 0, 0, 0, 0]
    # print(qs_df)
    # print(file_path)
    calc_score(df=qs_df, file=file_path)
    
    # # df.loc[df["#Measure"] == "Total Sequences", "Value"])

    chart_coloring = getPASScolors(summary_df)
    io.renderers.default="svg"
    makeOtherRadarChart(df=qs_df, file=file_path, colors=chart_coloring)
