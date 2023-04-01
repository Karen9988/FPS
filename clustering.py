import os, sys
import numpy as np
import pandas as pd
import csv

if len(sys.argv) < 4:
    print("Usage: python match.py {cadidate_path} {FPS_POS_path} {FPS_NEG_path} ")
    sys.exit(-1)

"""
Input data in CSV file
"""
cadidate_file = sys.argv[1]
pos_result_file = sys.argv[2]
neg_result_file = sys.argv[3]

# 1. read all file content
## 1.1 read cadidate
cadidate = pd.read_csv(cadidate_file)
# print(cadidate)
# cadidate may look like:
#       ID         Rt    Mz_POS  Size
# 0    330   3.653100  417.1365    61
# 1    345   3.745017  465.1032   104
# ......

## 1.2 read FPS_POS
fps_pos = pd.read_csv(pos_result_file)

## 1.3 read FPS_NEG
fps_neg = pd.read_csv(neg_result_file)

# 2. now we get the MZ_POS value of cadidates
candidate_mz_pos = cadidate["Mz_POS"]
# print(candidate_mz_pos)

may_same_collector_pos = []
# 3. Then we match cadidates from FPS_POS with mz_pos value
for candidate_idx, mz_pos in candidate_mz_pos.iteritems():
    # calc mz_pos resonable value range
    threshold = 0.05   # may passed by argument
    lower_limit = mz_pos - threshold
    upper_bound = mz_pos + threshold
    # print("lower_limit:%f, upper_bound: %f" % (lower_limit, upper_bound))
    may_same = []
    # get FPS_POS mz_pos data
    fps_pos_values = fps_pos[["Precursor m/z", "RT (min)", "Height"]]
    for idx, pos_row in fps_pos_values.iterrows():
        mz_value = pos_row["Precursor m/z"]
        if mz_value >= lower_limit and mz_value <= upper_bound:
            may_same.append(pos_row)
        # may_same_collector.append(may_same)
    def sort_by_rt(row):
        return row["RT (min)"]
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
    threshold = 0.05   # may passed by argument
    lower_limit = mz_neg - threshold
    upper_bound = mz_neg + threshold
    # print("lower_limit:%f, upper_bound: %f" % (lower_limit, upper_bound))
    may_same = []
    # get FPS_NEG mz_neg data
    fps_neg_values = fps_neg[["Precursor m/z", "RT (min)", "Height"]]
    for idx, neg_row in fps_neg_values.iterrows():
        mz_value = neg_row["Precursor m/z"]
        if mz_value >= lower_limit and mz_value <= upper_bound:
            may_same.append(neg_row)
        # may_same_collector.append(may_same)
    def sort_by_rt(row):
        return row["RT (min)"]
    may_same.sort(key=sort_by_rt)
    may_same_collector_neg.append(may_same)
# print(may_same_collector_neg)
# print("//////////////////////////")
# 5. compare candidate_pos and candidate_neg by RT,
#    select what we expected
matched_pos_and_neg = []
rt_threshold = 0.5
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
        rt_diff = abs(pos_candidates[pos_index]["RT (min)"] - neg_candidates[neg_index]["RT (min)"])
        if rt_diff < rt_threshold:
            # pos and neg are matched
            may_match.append((pos_candidates[pos_index], neg_candidates[neg_index]))
            pos_index += 1
            neg_index += 1
        else:
            if pos_candidates[pos_index]["RT (min)"] < neg_candidates[neg_index]["RT (min)"]:
                pos_index += 1
            else:
                neg_index += 1
    matched_pos_and_neg.append(may_match)
# print(matched_pos_and_neg)
# print("###############################")

# 6. calc pos_height/neg_height
pos_height_divide_neg_height = []
for pos_and_neg_list in matched_pos_and_neg:
    result_list = []
    for pos_and_neg in pos_and_neg_list:
        ratio = pos_and_neg[0]["Height"] / pos_and_neg[1]["Height"]
        if ratio < 1:
            ratio = -1 / ratio
        result_list.append(ratio)
    pos_height_divide_neg_height.append(result_list)
# print(pos_height_divide_neg_height)

# 7. cluster into 'Anthocyanin' and 'Flavonol glycoside'
## 'Anthocyanin': > 4.5
## 'Flavonol glycoside': < -1
candidate_idx = 0
header = ['candidate_index', 'Ratio', 'Cluster', 'pos_id', 'pos_mz', 'pos_rt',
          'pos_height', 'neg_id', 'neg_mz', 'neg_rt', 'neg_height']
data = []
for value_list in pos_height_divide_neg_height:
    # print("Matched Results for Candidate#%d" % candidate_idx)
    pos_and_neg_list = matched_pos_and_neg[candidate_idx]
    if len(value_list) == 0:
        print("No Matched Cands found!")
    value_idx = 0
    for value in value_list:
        pos_item = pos_and_neg_list[value_idx][0]
        neg_item = pos_and_neg_list[value_idx][1]
        value_idx += 1
        cluster = 'Unknown'
        if value > 4.5:
            cluster = 'Anthocyanin'
        elif value < -1:
            cluster = 'Flavonol glycoside'
        data.append([candidate_idx, value, cluster, pos_item.name, pos_item["Precursor m/z"], pos_item["RT (min)"],
                     pos_item["Height"], neg_item.name, neg_item["Precursor m/z"], neg_item["RT (min)"], neg_item["Height"]])
        # print("Ratio:%f, cluster:%s, pos_info:{\"id\":%s, \"mz\":%f, \"rt\":%f, \"height\":%f},"
        #     "neg_info:{\"id\":%s, \"mz\":%f, \"rt\":%f, \"height\":%f}" % 
        #     (value, cluster, pos_item.name, pos_item["Precursor m/z"], pos_item["RT (min)"],
        #     pos_item["Height"], neg_item.name, neg_item["Precursor m/z"], neg_item["RT (min)"],
        #     neg_item["Height"]))
    candidate_idx += 1
    print("")
    
# 8. Write data to CSV file
with open('result.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(data)
