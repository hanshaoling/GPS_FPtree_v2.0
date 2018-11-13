import math
import datetime
from networkx import all_neighbors, ancestors, descendants
import dill
load_start = datetime.datetime.now()
dill_file = 'base.pkl'
dill.load_session(dill_file)
load_end = datetime.datetime.now()
print('Load time: ', load_end - load_start)


threshold = input('The minimal SUPPORT is set to: ')
if threshold == '':
    print('Default value: 10')
    threshold = 10
else:
    threshold = int(threshold)
start1 = datetime.datetime.now()
FPTREE, HEADERTABLE = create_fp_tree(init_dict, threshold)

# Find the parent of a given node and output a complete path.
# This function will be called and explained in the next function.
# prefix_path will be assigned as an empty list.


def complete_path(node, prefix_path):
    # node is an object in class node and prefix_path is a list.
    if node.parent is not None:
        prefix_path.append(node.name)
        complete_path(node.parent, prefix_path)


def find_complete_path(phenotype, header_table):
    node = header_table[phenotype][1]
    # A node class object, as the lowest level of this path, corresponding a certain phenotype.
    complete_set = {}
    # a new dict to store the results: a set of all possible group of phenotypes starting from a certain phenotype.
    # (searching in one direction)
    while node is not None:
        prefix_path = []
        complete_path(node, prefix_path)
        # results are stored in list 'prefix_path'
        if len(prefix_path) > 1:
            complete_set[frozenset(prefix_path[1:])] = node.count
        node = node.nodeLink
        # continue to the next node(but with the same phenotype)
    return complete_set
    # a dict containing key : value pair -- frozenset of prefix of a phenotype :  count of this prefix.

# Capsule these functions.


def mine_fp_tree(FPTREE, header_table, min_sup = 10, freq_item_list = [], prefix = set([])):
    ordered = sorted(header_table.items(), key=lambda x: x[0])
    ordered_header_table = []
    for each in ordered:
        ordered_header_table.append(each[0])

    for each in ordered_header_table:
        new_freq_set = prefix.copy()
        new_freq_set.add(each)
        freq_item_list.append(new_freq_set)
        condition_base = find_complete_path(each, header_table)
        condition_tree, condition_head = create_fp_tree(condition_base, min_sup)
        if condition_head is not None:
            mine_fp_tree(condition_tree, condition_head, min_sup, freq_item_list, new_freq_set)
    return freq_item_list


final = mine_fp_tree(FPTREE, HEADERTABLE, threshold)
sorted_list = sorted(final, key=lambda x: len(x), reverse=True)  # Sort the list by length.


def get_subset(set):  # Get all REAL subset of a listed-freq_set. Return is a list of list. [[...], [...], [...], ...]
    N = len(set)
    subset = []
    for i in range(2 ** N):
        combo = []
        for j in range(N):
            if ( i >> j ) % 2 == 1:
                combo.append(set[j])
        if combo != [] and len(combo) != N:
            subset.append(combo)
    return subset


memo = {}  # Counting for subset.


def count_subset(subset):  # subset is a frozenset derived from list.
    ls = list(init_dict.keys())
    count = 0
    if subset in memo:
        return memo[subset]
    for each in ls:
        if subset.issubset(each):
            count += 1
    memo[subset] = count
    return count


def get_association(sorted_list, min_confidence=0.5):
    rules = {}
    for freq_set in sorted_list:
        freq_listed = list(freq_set)
        freq_subset = get_subset(freq_listed)
        for each in freq_subset:
            a = frozenset(each)
            count = count_subset(a)
            confidence = count_subset(frozenset(freq_set)) / count
            if confidence >= min_confidence:
                good_rule = (a, frozenset(freq_set) - a)
                rules[good_rule] = confidence
    return rules


sorted_list_origin = []  # Convert the output from index to origin phenotype name.
for i in range(len(sorted_list)):
    sorted_list_origin.append(list(sorted_list[i]))
    for j in range(len(sorted_list[i])):
        sorted_list_origin[i][j] = term_name_dict[int(sorted_list_origin[i][j])]
sorted_list_origin = sorted(sorted_list_origin, key=lambda x: len(x), reverse=True)


confidence_threshold = input('The minimal CONFIDENCE is set to: ')
if confidence_threshold == '':
    print('Default value: 0.5')
    confidence_threshold = 0.5
else:
    confidence_threshold = float(confidence_threshold)
rules = get_association(sorted_list, confidence_threshold)



def lift_filter(rule_from, rule_to):
    n_total = len(set(disease_origin))
    lift = (memo[frozenset(list(rule_from) + list(rule_to))] * n_total) / (memo[rule_from] * memo[rule_to])
    return lift


for rule, confidence in rules.items():
    rule_list = list(rule)
    lift = lift_filter(rule_list[0], rule_list[1])
    rules[rule] = (confidence, lift)


AnnoNum = {}
# Dict, phenoID : number of annotation
for each in convert.values():
    each = list(each)
    for every in each:
        if every not in AnnoNum:
            AnnoNum[every] = 1
        if every in AnnoNum:
            AnnoNum[every] += 1

end1 = datetime.datetime.now()
print('Run time: ', end1 - start1)


def get_ic(term):
    return math.log(AnnoNum[1]/AnnoNum[term])


def get_mica_ic(term1, term2):
    ancestor_set1 = ancestors(G, term1)
    ancestor_set1.add(1)
    ancestor_set2 = ancestors(G, term2)
    ancestor_set2.add(1)
    common_set = ancestor_set1 & ancestor_set2
    sorted_list = sorted(list(common_set), key=lambda x: get_ic(x), reverse=True)
    return get_ic(sorted_list[0])


def similarity(term1, term2):
    sim = 2 * get_mica_ic(term1, term2) / (get_ic(term1) + get_ic(term2))
    return sim
# value is in interval (0,1)


def set_similarity(set1, set2):
    ls = []
    set_sim = 0
    for A in set1:
        for B in set2:
            ls.append((A, B))
    for each in ls:
        set_sim += similarity(each[0], each[1]) ** 2
    return math.sqrt(set_sim / len(ls))


while True:
    set_input = input('phenotype input: ').split(sep = ',')
    start = datetime.datetime.now()
    for i in range(len(set_input)):
        set_input[i] = int(set_input[i])
    set_input = frozenset(set_input)
    set_input_name = set()
    for each in set_input:
        set_input_name.add(term_name_dict[each])
    rules_filter = {}
    for rule, tuple in rules.items():
        a = set_similarity(set_input, rule[0])
        # tuple: confidence, lift, memo: support of FROM
        rules_filter[rule] = tuple + (a, memo[rule[0]])
    rules_filter = sorted(rules_filter.items(), key=lambda x: (x[1][2], x[1][0]), reverse=True)

    rules_wanted = []
    for x in rules_filter[:5]:
        FROM_LIST = list(x[0][0])
        TO_LIST = list(x[0][1])
        for i in range(len(FROM_LIST)):
            FROM_LIST[i] = term_name_dict[FROM_LIST[i]]
        for j in range(len(TO_LIST)):
            TO_LIST[j] = term_name_dict[TO_LIST[j]]
        rules_wanted.append([FROM_LIST, TO_LIST, list(x[1])])

    '''for each_rule in rules_wanted:
        print('From: ', each_rule[0], '\n', 'To: ', each_rule[1], '\nConfidence: ', each_rule[2][0], '|Lift: ', each_rule[2][1], '|Similarity: ', each_rule[2][2], '|Support: ', each_rule[2][3], '\n')
'''
    waiting_list = []
    for each_rule in rules_wanted:
        waiting_list += each_rule[0]
        waiting_list += each_rule[1]
    print(set(waiting_list) - set_input_name)
    print(len(set(waiting_list) - set_input_name))
    end = datetime.datetime.now()
    print('Search time: ', end - start)
    print('###############################################################################################' * 3)
