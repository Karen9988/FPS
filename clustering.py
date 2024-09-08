# -*- coding: UTF-8 -*-
import os, sys
import numpy as np
import pandas as pd
import csv

if len(sys.argv) < 4:
    print("Usage: python clustering.py {candidate_path} {FPS_POS_path} {FPS_NEG_path} ")
    sys.exit(-1)

"""
Input data in CSV file
"""
candidate_file = sys.argv[1]
pos_result_file = sys.argv[2]
neg_result_file = sys.argv[3]

# 1. read all file content
## 1.1 read candidate
candidate = pd.read_csv(candidate_file)
# print(candidate)
# candidate may look like:
#       ID         Rt    Mz_POS  Size
# 0    330   3.653100  417.1365    61
# 1    345   3.745017  465.1032   104
# ......

## 1.2 read FPS_POS
fps_pos = pd.read_csv(pos_result_file)

## 1.3 read FPS_NEG
fps_neg = pd.read_csv(neg_result_file)

# 2. now we get the MZ_POS value of candidates
candidate_mz_pos = candidate["Mz_POS"]
# print(candidate_mz_pos)

may_same_collector_pos = []
# 3. Then we match candidates from FPS_POS with mz_pos value
for candidate_idx, mz_pos in candidate_mz_pos.iteritems():
    # calc mz_pos resonable value range
    tolerance = 0.01   # may passed by argument
    lower_bound = mz_pos - tolerance
    upper_bound = mz_pos + tolerance
    #print("lower_bound: %f, upper_bound: %f" % (lower_bound, upper_bound))
    may_same = []
    # get FPS_POS mz_pos data
    fps_pos_values = fps_pos[["Precursor m/z", "RT (min)", "Area"]]
    for idx, pos_row in fps_pos_values.iterrows():
        mz_value = pos_row["Precursor m/z"]
        if mz_value >= lower_bound and mz_value <= upper_bound:
            may_same.append((idx, pos_row, abs(mz_pos - mz_value)))
        # may_same_collector.append(may_same)
    def sort_by_rt(row):
        return row[1]["RT (min)"]
    may_same.sort(key=sort_by_rt)
    may_same_collector_pos.append(may_same)
# print(may_same_collector_pos)
# print("//////////////////////////")

# 4. then match candidates from FPS_NEG with mz_pos value
may_same_collector_neg = []
## 4.1 candidate_neg = candidate_pos - 2H
candidate_mz_neg = candidate_mz_pos - 2 * 1.007825
## 4.2 select candidate
for candidate_idx, mz_neg in candidate_mz_neg.iteritems():
    # calc mz_neg resonable value range
    tolerance = 0.01   # may passed by argument
    lower_bound = mz_neg - tolerance
    upper_bound = mz_neg + tolerance
    # print("lower_bound:%f, upper_bound: %f" % (lower_bound, upper_bound))
    may_same = []
    # get FPS_NEG mz_neg data
    fps_neg_values = fps_neg[["Precursor m/z", "RT (min)", "Area"]]
    for idx, neg_row in fps_neg_values.iterrows():
        mz_value = neg_row["Precursor m/z"]
        if mz_value >= lower_bound and mz_value <= upper_bound:
            may_same.append((idx, neg_row, abs(mz_value - mz_neg)))
        # may_same_collector.append(may_same)
    def sort_by_rt(row):
        return row[1]["RT (min)"]
    may_same.sort(key=sort_by_rt)
    may_same_collector_neg.append(may_same)
# print(may_same_collector_neg)
# print("//////////////////////////")
# 5. compare candidate_pos and candidate_neg by RT,
#    select what we expected
matched_pos_and_neg = []
rt_tolerance = 0.02
for i in range(0, len(may_same_collector_pos)):
    pos_candidates = may_same_collector_pos[i]
    pos_index = 0
    pos_end = len(pos_candidates)
    neg_candidates = may_same_collector_neg[i]
    neg_index = 0
    neg_end = len(neg_candidates)
    may_match = []
    while (pos_index < pos_end) and (neg_index < neg_end):
        # calc rt-diff
        rt_diff = abs(pos_candidates[pos_index][1]["RT (min)"] - neg_candidates[neg_index][1]["RT (min)"])
        if rt_diff < rt_tolerance:
            # pos and neg are matched
            already_found_pair = False
            for matched_array in matched_pos_and_neg:
                for matched in matched_array:
                    if pos_candidates[pos_index][0] == matched[0][0] and neg_candidates[neg_index][0] == matched[1][0]:
                        already_found_pair = True
                        break
                if already_found_pair:
                    break
            if not already_found_pair:
                may_match.append((pos_candidates[pos_index], neg_candidates[neg_index], rt_diff))
            pos_index += 1
            neg_index += 1
        else:
            if pos_candidates[pos_index][1]["RT (min)"] < neg_candidates[neg_index][1]["RT (min)"]:
                pos_index += 1
            else:
                neg_index += 1
    matched_pos_and_neg.append(may_match)
# print(matched_pos_and_neg)
# print("###############################")

# 6. calc pos_area/neg_area
pos_area_divide_neg_area = []
for pos_and_neg_list in matched_pos_and_neg:
    result_list = []
    for pos_and_neg in pos_and_neg_list:
        ratio = pos_and_neg[0][1]["Area"] / pos_and_neg[1][1]["Area"]
        rt_diff = pos_and_neg[2]
        if ratio < 1:
            ratio = -1 / ratio
        result_list.append((ratio, rt_diff))
    pos_area_divide_neg_area.append(result_list)
# print(pos_area_divide_neg_area)

# 7. cluster into 'Anthocyanin' and 'Flavonol glycoside'
## 'Anthocyanin': > 4.5
## 'Flavonol glycoside': < -1
candidate_idx = 0
header = ['Candidate_index', 'POS_RT', 'POS_Precursor', 'POS_difference', 'POS_Area',
          'NEG_RT', 'RT_difference', 'NEG_Precursor', 'NEG_difference', 'NEG_Area', 'Ratio', 'Cluster']
data = []
for value_list in pos_area_divide_neg_area:
    # print("Matched Results for Candidate#%d" % candidate_idx)
    pos_and_neg_list = matched_pos_and_neg[candidate_idx]
    if len(value_list) == 0:
        print("No Matched Candidates found!")
    value_idx = 0
    for ratio_and_rtdiff in value_list:
        pos_item = pos_and_neg_list[value_idx][0][1]
        neg_item = pos_and_neg_list[value_idx][1][1]
        cluster = 'Unknown'
        value = ratio_and_rtdiff[0]
        if value > 5:
            cluster = 'A'
        elif value < -1:
            cluster = 'F'
        data.append([candidate_idx, pos_item["RT (min)"], pos_item["Precursor m/z"], pos_and_neg_list[value_idx][0][2],
                     pos_item["Area"], neg_item["RT (min)"], ratio_and_rtdiff[1], neg_item["Precursor m/z"],
                     pos_and_neg_list[value_idx][1][2], neg_item["Area"], value, cluster])
        value_idx += 1
    candidate_idx += 1
    # print("")

# 8. Write data to CSV file
with open('result.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(data)
