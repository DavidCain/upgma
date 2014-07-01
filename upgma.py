#!/usr/bin/env python

"""
Builds an UPGMA tree- a taxonomy tree built on averaged taxa distances.
"""

import itertools


class Node(object):
    def __init__(self, dist_function, left=None, right=None):
        """ Initialize a taxon or OTU.

        If Node is to be a single taxon, set left, leave right==None.
        If Node is an OTU, set left and right to other Node objects.

        dist_function: function to calculate the distance between
        two Residue objects. Must take two parameters and return a
        numeric distance.
        """
        self.dist_function = dist_function
        self.left = left
        self.right = right
        self.nn = None

    def get_nn(self):
        """ Return the nearest neighbor found by the last operation. """
        return self.nn

    def get_left(self):
        return self.left

    def get_right(self):
        return self.right

    def __str__(self):
        return "-".join(str(node) for node in self.__iter__())

    def __iter__(self):
        """ Yield all taxa in the Node and its children.

        Will not yield OTU's, only actual taxa.

        The iterator is useful for averaging distances, and returning
        all members of a branch.
        """
        for child in ([self.left, self.right]):
            if child is None:
                pass
            # TODO: Not a Pythonic approach to solving the problem
            elif isinstance(child, Node):
                for item in child:
                    yield item
            else:
                yield child

    def __len__(self):
        """ Returns how many taxa are contained within the Node, (not
        counting OTU's, just taxa)
        """
        return sum(1 for taxa in self.__iter__())

    def get_dist(self, node2):
        """ Get the arithmetic mean distance between the node and another.

        For finding the distance between two childless taxa nodes, uses
        the distance function defined in __init__ (then averages these
        distances ).
        """
        res_pairs = itertools.product(self, node2)
        dist = sum(self.dist_function(*pair) for pair in res_pairs)
        return dist / (len(self) + len(node2))

    def update_distances(self, clusters):
        """ Update distances, find and store the nearest neighbor.

        Returns the nearest neighbor
        TODO: store distances in heap? (see paper)
        """
        # Find NN, and distance to it
        self.nn = None
        self.dist_to_NN = None
        for node in clusters:
            if node is self:
                continue
            dist_to_node = self.get_dist(node)
            no_NN = (self.dist_to_NN is None)
            if no_NN or (dist_to_node < self.dist_to_NN):
                self.dist_to_NN = dist_to_node
                self.nn = node
        return self.nn


class UPGMA(object):
    def __init__(self, taxa, dist_function):
        """ Build an efficient UPGMA tree.

        Build a UPGMA tree using the O(n^2) approach detailed in: Gronau
        I, Moran S (2007) Optimal implementations of UPGMA and other
        common clustering algorithms.

        Nodes are clustered according to their proximity, determined by
        the distance function (function must take two parameters and
        return a numeric distance).
        """
        clusters = set([Node(dist_function, taxon) for taxon in taxa])
        self.dist_function = dist_function
        self.build_tree(clusters)

    def get_tree(self):
        """ Return the UPGMA tree.

        The "tree" is just a pointer to the root Node object, which
        contains all the other Nodes as its children.
        """
        return self.tree

    def get_largest_branch(self):
        """ Return the branch with more final taxa (not counting OTU's). """
        if self.tree.get_right() is None:
            #print "No right branch, returning left branch."
            return self.tree.get_left()

        elif self.tree.get_left() is None:
            #print "No left branch, returning right branch."
            return self.tree.get_right()

        elif len(self.tree.get_right()) > len(self.tree.get_left()):
            return self.tree.get_right()

        else:
            return self.tree.get_left()

    def build_tree(self, clusters):
        """ Build the UPGMA tree.

        Goes through clusters, grouping the two closest together into an
        OTU. Grouping continues until all clusters belong to an OTU.
        """
        for node in clusters:
            node.update_distances(clusters)

        for i in range(len(clusters) - 1):
            #print 18*"-"
            # Find globally closest pair
            c1, c2 = self._get_closest_pair(clusters)

            # Join the clusters into an OTU
            clusters.remove(c1)
            clusters.remove(c2)
            new_node = self._create_OTU(c1, c2)
            clusters.add(new_node)

            #print "Finding nearest neighbor for new cluster"
            new_node.update_distances(clusters)

            # Recalculate NN's for each cluster previously having one of the
            # joined clusters as its NN
            for node in clusters:
                if (node.get_nn() == c1) or (node.get_nn() == c2):
                    node.update_distances(clusters)

        # At this point "clusters" contains just the root of the tree
        assert len(clusters) == 1, "UPGMA tree construction failed"
        self.tree = clusters.pop()

    def _create_OTU(self, c1, c2):
        """ Join two clusters into an OTU. """
        return Node(self.dist_function, c1, c2)

    def _get_closest_pair(self, clusters):
        """ Return the two globally closest clusters of a collection. """
        min_dist = None
        for node in clusters:
            if (min_dist is None) or (node.dist_to_NN < min_dist):
                min_dist = node.dist_to_NN
                c1 = node
                c2 = node.get_nn()
                assert node.get_nn() is not None
        #print "Clusters {} and {} closest".format(c1, c2)
        return (c1, c2)


if __name__ == '__main__':
    def simple_diff(o1, o2):
        """ Return a simple absolute distance between the objects. """
        return abs(o1 - o2)

    from random import randint
    random_integers = [randint(0, 200) for i in range(8)]
    print "Initial clusters: {}".format(random_integers)
    print UPGMA(random_integers, simple_diff).get_tree()
