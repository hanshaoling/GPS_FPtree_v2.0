import dill
import re
import networkx as nx
import datetime
start = datetime.datetime.now()
# deal with HPO terms
filehandle = open('hpo1009.txt')
fileload = filehandle.read()
txtsplit = fileload.split('\n\n')
# txtsplit is a list of chunks.
terms = txtsplit[1:]

term_id_list = []  # including obsolete terms
term_id_list_inuse = []  # only with terms that have parent
term_name_list = []
term_name_dict = {}  # matching HPO term id and name
parent_offspring = set()  # record relationship of parent-offspring, tuple(parent's id, offspring's id)

for term in terms:
    id = int(re.search('id: HP:(\d+)', term).group(1))
    name = re.search('\nname: (.+)\n', term).group(1)
    parents = re.findall('is_a: HP:(\d+)', term)
    if len(parents) != 0:
        term_id_list_inuse.append(id)
        for i in range(len(parents)):
            parents[i] = int(parents[i])
    term_id_list.append(id)
    term_name_list.append(name)
    if len(parents) != 0:
        for i in range(len(parents)):
            parent_offspring.add((parents[i], id))

for i in range(len(term_id_list)):
    term_name_dict[term_id_list[i]] = term_name_list[i]

list_PO = list(parent_offspring)  # list of pairing relationship
G = nx.DiGraph()
G.add_edges_from(list_PO)  # drawing DAG


def inflate(ls):  # to inflate the phenotype set of each disease, ls is a list
    for each in ls:
        ls.extend(list(nx.ancestors(G, each)))
    return set(ls)


# deal with annotation
filename = 'FPTREE.txt'  # the file has been cleaned up, replacing the obsolete terms
txt = open(filename)
data = txt.readlines()
data = data[1:]
# Open and load the file. 'data' is a list.
disease_origin = []
phenotype_origin = []

for x in data:
    disease_i = x.split('\t')[0].strip()
    phenotype_i = x.split('\t')[1].strip()
    disease_origin.append(int(disease_i))
    phenotype_origin.append(int(phenotype_i))

convert = {}  # matching disease and related HPO terms (inflated), disease's id : [list of HPO terms]
length = len(disease_origin)
for i in range(length):
    if disease_origin[i] not in convert:
        convert[disease_origin[i]] = [phenotype_origin[i]]
    else:
        convert[disease_origin[i]].append(phenotype_origin[i])
for each_disease in list(convert.keys()):
    convert[each_disease] = inflate(convert[each_disease])  # directly change the value in dict(convert)
phenotype_set = list(convert.values())


class Node:  # building fp-tree
    def __init__(self, name, count, parent):
        # 5 attributes.
        self.name = name
        self.count = count
        self.parent = parent
        self.nodeLink = None
        self.children = {}
        # The nodes next level to self, key : value is name : node.

    def increase(self, count):
        self.count += count
        # counter of the node


def create_init(phenotype_set):  # This function turn the list of lists into a dict counting each group of phenotypes.
    init_dict = {}  # Dict[FrozenSet[Any], int], in this case the int should all be 1 (1 stands for each disease)
    for event in phenotype_set:
        key = frozenset(event)
        if key in init_dict:
            init_dict[key] += 1
        else:
            init_dict[key] = 1
    return init_dict


def update_header(node_to_test, target_node):  # node_to_test belongs to class Node
    while node_to_test.nodeLink is not None:
        node_to_test = node_to_test.nodeLink
    node_to_test.nodeLink = target_node


def update_fp_tree(event, in_tree, header_table, count):  # This function will be needed in the next function.
    # event is a list of phenotypes of one disease,
    # in_tree is a destination node, in class node.
    # header_table is just header_table.
    # count is the occurring number of this event. (see create_init())
    if event[0] in in_tree.children:
        # check if the first item of event is already a child of inTree.
        in_tree.children[event[0]].increase(count)
        # children is also in class node.
    else:
        # create new branch.
        in_tree.children[event[0]] = Node(event[0], count, in_tree)
        if header_table[event[0]][1] is None:
            header_table[event[0]][1] = in_tree.children[event[0]]
        else:
            update_header(header_table[event[0]][1], in_tree.children[event[0]])
    if len(event) > 1:
        update_fp_tree(event[1:], in_tree.children[event[0]], header_table, count)


def create_fp_tree(init_dict, min_sup):  # The main executor of creating fptree.
    header_table = {}
    for event in init_dict:
        for phenotype in event:
            header_table[phenotype] = header_table.get(phenotype, 0) + init_dict[event]
            # Counting for every single phenotype.
    for k in list(header_table.keys()):
        if header_table[k] < min_sup:
            del (header_table[k])
            # Delete the phenotypes with occurring number < min_sup.
    freq_pheno_set = set(header_table.keys())
    # Make a set of the 'freq' phenotypes.
    if len(freq_pheno_set) == 0:
        # min_sup is too large.
        return None, None
    for k in header_table:
        header_table[k] = [header_table[k], None]
        # Change the format: element [count, node].

    root_tree = Node('Null Set', 1, None)
    # Make a new ROOT node, with name 'Null Set', count '1', and no parent.
    for event, count in init_dict.items():
        # init_dict：[element(group of phenotypes), count]
        select_phenotype = {}
        for phenotype in event:
            if phenotype in freq_pheno_set:
                select_phenotype[phenotype] = header_table[phenotype][0]
                # Select the phenotypes with support > minSup.
                # In select_phenotype, phenotype : count
        if len(select_phenotype) > 0:
            # Order the phenotypes by counting number.
            ordered_list = sorted(select_phenotype.items(), key=lambda x: x[1], reverse=True)
            ordered_phenotype = []
            for each in ordered_list:
                ordered_phenotype.append(each[0])
            # Return a list of phenotype, selected and ordered.
            update_fp_tree(ordered_phenotype, root_tree, header_table, count)
    return root_tree, header_table


init_dict = create_init(phenotype_set)

annotation_number = {}  # Dict, HPO term's ID : number of annotation. for calculation of Information Content
for each in convert.values():
    each = list(each)
    for every in each:
        if every not in annotation_number:
            annotation_number[every] = 1
        if every in annotation_number:
            annotation_number[every] += 1

dill_file = 'base.pkl'
dill.dump_session(dill_file)
end = datetime.datetime.now()
print('Run time: ', end - start) 