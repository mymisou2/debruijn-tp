#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""
from operator import itemgetter
import argparse
import os
import sys
import random
random.seed(9001)
import statistics
import matplotlib as plt
import networkx as nx

__author__ = "Myriam Mendli"
__copyright__ = "EISTI"
__credits__ = ["Myriam Mendli"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Myriam Mendli"
__email__ = "mendlimyri@eisti.eu"
__status__ = "Rendu"

def isfile(path):
    """
        Check if path is an existing file.
        :Parameters:
        path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """
    Retrieves the arguments of the program.
    Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    """
        Genreateur de séquences
        Input:
            fastq_file: fichier fastq
        Return
            un générateur de séquences
    """
    with open(fastq_file, "rt") as monf:
        for line in monf:
            yield next(monf).replace("\n", "")
            next(monf)
            next(monf)

def cut_kmer(read, kmer_size):
    """
        Générateur de k-mer
        Input:
            read: une séquence
            kmer_size : une taille de k-mer
        Output:
            un générateur de k-mer
    """
    i = 0
    tmp = i + kmer_size
    while tmp < len(read) + 1:
        kmer = read[i:i+kmer_size]
        yield kmer
        i = i+1
        tmp = i + kmer_size


def build_kmer_dict(fastq_file, kmer_size):
    """
        Générateur de k-mer
        Input:
            read: une séquence
            kmer_size : une taille de k-mer
        Output:
            un générateur de k-mer
    """
    if isfile(fastq_file):
        dico = dict()
        for generator_file in read_fastq(fastq_file):
            k_mers = cut_kmer(generator_file, kmer_size)
            for k_mer in k_mers:
                if k_mer in dico:
                    dico[k_mer] += 1
                else:
                    dico[k_mer] = 1
        return dico
    return None

############################################################

def build_graph(kmer_dict):
    """ Build a graph using a k-mer dictionnary """
    graph = nx.DiGraph()
    for key in kmer_dict.keys():
        graph.add_edge(key[:-1], key[1:], weight=kmer_dict[key])
    return graph

def get_starting_nodes(graph):
    """ prend en entrée un graphe et retourne une liste de noeuds d’entrée
     soit de degré d'entré 0 """
    return [n for n, d in graph.in_degree() if d == 0]

def get_sink_nodes(graph):
    """ prend en entrée un graphe et retourne une liste de noeuds
     de sortie soit de degré de sorti 1 ou 0 """
    return [x for x in graph.nodes() if graph.out_degree(x) == 0 and graph.in_degree(x) == 1]

def get_contigs(graph, starting_nodes, ending_nodes):
    """ prend un graphe, une liste de noeuds d’entrée et une liste de
     sortie et retourne une liste de tuple(contig, taille du contig) """
    # Calcul distance
    distances = dict(nx.all_pairs_shortest_path_length(graph))
    tuple = []
    for start_node in starting_nodes:
        for end_node in ending_nodes:
            path = ""
            # Foreach distance from the begininng to the end
            for path_end2end in nx.all_simple_paths(graph, start_node, end_node):
                # Take only the fist element of the list
                concatention = [value for key, value in enumerate(path_end2end, 1) if key%2 != 0]
                # Create a path with all the element
                for element in concatention:
                    path += element
            if path != "":
                # Add the tuple contig  with distance
                tuple.append((path, distances[start_node][end_node]+2))
    return tuple

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(contigs_list, output_file):
    """  qui prend une liste de tuple (contig, taille du contig) et un nom de
     fichier de sortie et écrit un fichier de sortie contenant les
     contigs selon le format fasta (retour chariot tous les 80 caractères
    à l’aide de la fonction fill: """
    with open(output_file, "w+") as file:
        for i, _ in enumerate(contigs_list):
            file.write(">contig_" + str(i) +" len=" + str(contigs_list[i][1]) + "\n")
            fill_ = fill(contigs_list[i][0])
            file.write(fill_+"\n")

############################################################

def std(data):
    """ prend une liste de valeur, qui retourne l’écart type """
    return statistics.stdev(data)

def path_average_weight(graph, path):
    """ prend un graphe et un chemin et qui retourne un poids moyen """
    moyenne = 0.0
    nb_arcs = 0
    for indice_node, _ in enumerate(path):
        if indice_node != len(path) - 1:
            try:
                moyenne += graph[path[indice_node]][path[indice_node + 1]]['weight']
                nb_arcs += 1
            except:
                pass
    try:
        moyenne = moyenne/nb_arcs
    except:
        return 0
    return moyenne

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """ prend un graphe et une liste de chemin, la variable booléenne """
    # Pour chaque chemin
    for path in path_list:
        remove = list()
        # On considère le premier element et le dernier element du chemin
        starting_node = path[0]
        ending_node = path[-1]

        for node in path:
            # Pour chaque noeud on test s'il est un noeud début ou fin
            # pour savoir si on peut le supprimer
            if(delete_entry_node is False and node == starting_node):
                pass
            else:
                if(delete_sink_node is False and node == ending_node):
                    pass
                else:
                    # On update la liste de noeuds a supprimer
                    remove.append(node)
        # On supprime les noeuds a supprmier
        graph.remove_nodes_from(remove)

    return graph

def index_weight(weight_avg_list, path_length):
    """ Return the index weight """
    # Recupère les index des chemins ayant les poinds le poid le plus lourd
    best_weight_index = [i for i, j in enumerate(weight_avg_list) if j == max(weight_avg_list)]
    index_weight_keep = best_weight_index[0]

    # S'il y en aplusieurs on selectionne par la taille du chemin
    if len(best_weight_index) > 1:
        exist_longest = True
        size_lg = path_length[index_weight_keep]
        # Select le plus long
        # On parcourt les index des elements ayant les poids les plus lourds
        for i in best_weight_index[1:]:
            # Si 2 elements ont la meme taille
            if size_lg == path_length[i]:
                exist_longest = False
            else:
                # Si on trouve un max
                if path_length[i] > size_lg:
                    size_lg = path_length[i]
                    index_weight_keep = i
                    exist_longest = True
        # Si 2 chemins ont la meme taille
        if exist_longest is False:
            # Select random
            index_weight_keep = random.choice(best_weight_index)
    return index_weight_keep

def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    """ t retourne un graphe nettoyé des chemins indésirables """
    # On recupère l'index du chemin que l'on  va conserver
    index_weight_keep = index_weight(weight_avg_list, path_length)

    # Pour chque chemin de la liste
    for path in path_list:
        remove = list()
        # Si le chemin n'est pas le chemin que l'on souhaite conservé
        if path != path_list[index_weight_keep]:
            starting_node = path[0]
            ending_node = path[-1]
            # On efface chaque noeud qui ne sont pas des noeuds d'entrés ou de sorti
            for node in path:
                if(delete_entry_node is False and node == starting_node):
                    pass
                else:
                    if(delete_sink_node is False and node == ending_node):
                        pass
                    else:
                        remove.append(node)
        graph.remove_nodes_from(remove)
    return graph


#################################################################

def solve_bubble(graph, ancestor_node, descendant_node):
    """ Retourne un graph nettoyé de la bulle se trouvant entre ces deux nœuds """
    path_list = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    path_length = list()
    weight_avg_list = list()
    # Pour chaque chemin
    for path in path_list:
        # On calcul sa taille
        path_length.append(len(path))
        # On calucle son poids moyen
        weight_avg_list.append(path_average_weight(graph, path))
    # On renvoie le graph avec le meilleyur chemin
    return select_best_path(graph, path_list, path_length, weight_avg_list)

def simplify_bubbles(graph):
    """ Retourne un graphe sans bulle """
    nodes = list(graph.nodes)
    # POur chaque de noeud du graph
    for node in nodes:
        # Si le noeud n'a pas été supprimé du graph
        if node in graph.nodes:
            # On calcul ces predecessuers
            pred = list(graph.predecessors(node))
            # S'il eb a plusieurs
            if len(pred) > 1:
                for i, _ in enumerate(pred):
                    for j in range(i+1, len(pred)):
                        # On calcule l'ancetre commun le plus proche entre 2 noeuds
                        try:
                            ancestor = nx.algorithms.lowest_common_ancestor(graph, pred[i], pred[j])
                            # On recupère un graph nettoyé de la bulle
                            graph = solve_bubble(graph, ancestor, node)
                        except:
                            pass

    # On remet les noeuds supprimés
    [graph.add_node(node) for node in nodes if node not in graph.nodes]
    return graph

#################################################################


def solve_entry_tips(graph, starting_nodes):
    """ retourne graphe sans chemin d’entrée indésirable """
    tmp_node = "tmpnode"
    # Add a temporary node at the end to create a bull
    graph.add_node(tmp_node)
    # Add edges between the temporary node and starting nodes
    [graph.add_edge(tmp_node, node, weight=0) for node in starting_nodes]
    graph = simplify_bubbles(graph)
    # remove the temporary node
    graph.remove_node(tmp_node)
    return graph

def solve_out_tips(graph, ending_nodes):
    """ retourne graphe sans chemine testSortie indésirable """
    tmp_node = "tmpnode"
    # Add a temporary node at the end to create a bull
    graph.add_node(tmp_node)
    # Add edges between the  endings nodes and temporary node
    [graph.add_edge(node, tmp_node, weight=0) for node in ending_nodes]
    graph = simplify_bubbles(graph)
    # remove the temporary node
    graph.remove_node(tmp_node)

    return graph

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    kmer_dict = build_kmer_dict("/home/eisti/debruijn-tp/data/eva71_hundred_reads.fq", 30)
    graph = build_graph(kmer_dict)
    graph = simplify_bubbles(graph)
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    graph = solve_entry_tips(graph, starting_nodes)
    graph = solve_out_tips(graph, ending_nodes)
    contigs_list = get_contigs(graph, starting_nodes, ending_nodes)
    save_contigs(contigs_list, "/home/eisti/debruijn-tp/data/testSortie.fna")
    # Get arguments
    # args = get_arguments()

if __name__ == '__main__':
    main()
