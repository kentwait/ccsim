class PairwiseNeighborhood(object):
    """Simulate neighborhood interactions instead of explicitly simulating space.
    
    For example, in a linear neighborhood with 5 blocks with boundary adjustment,
    -> 0 <-> 1 <-> 2 <-> 3 <-> 4 <-
    There are 5 neighborhood interactions, (4,0), (0,1), (1,2), (2,3), (3,4)
    Instead of picking a block and then picking a neighbor, pick a block-block interaction instead.
    Thus the simulation becomes topology independent
    """
    
    def __init__(self, neighborhood:dict={}, node_metadata:dict={}, interaction_metadata:dict={}):
        self.neighborhood = neighborhood  # dict node_id_a: set of node_id_b
        self.node_metadata = node_metadata  # dict node_id: metadata
        self.interaction_metadata = interaction_metadata  # dict (node_id_a, node_id_b): metadata

    def add_node(self, nid, meta=None):
        if nid in self.neighborhood.keys():
            raise Exception(f'Node {nid} already exists')
        self.neighborhood[nid] = set([])
        if meta:
            self.node_metadata[nid] = meta

    def add_node_metadata(self, nid, meta):
        if nid not in self.neighborhood.keys():
            raise Exception(f'Node {nid} does not exist')
        self.node_metadata[nid] = meta

    def add_interaction(self, nid1, nid2, meta=None):
        # Check if nid1 and nid2 exist
        if nid1 not in self.neighborhood.keys():
            self.add_node(nid1)
        if nid2 not in self.neighborhood.keys():
            self.add_node(nid2)
        # Interactions are bidirectional
        self.neighborhood[nid1].add(nid2)
        self.neighborhood[nid2].add(nid1)
        if meta:
            self.interaction_metadata[(nid1, nid2)] = meta
            self.interaction_metadata[(nid2, nid1)] = self.interaction_metadata[(nid1, nid2)]

    def add_interaction_metadata(self, nid1, nid2, meta):
        # Check if nid1 and nid2 exist
        if nid1 not in self.neighborhood.keys():
            raise KeyError(f'Invalid node id: {nid1}')
        if nid2 not in self.neighborhood.keys():
            raise KeyError(f'Invalid node id: {nid2}')
        self.interaction_metadata[(nid1, nid2)] = meta
        self.interaction_metadata[(nid2, nid1)] = self.interaction_metadata[(nid1, nid2)]

    def remove_interaction(self, nid1, nid2):
        # Check if nid1 and nid2 exist
        if nid1 not in self.neighborhood.keys():
            raise KeyError(f'Invalid node id: {nid1}')
        if nid2 not in self.neighborhood.keys():
            raise KeyError(f'Invalid node id: {nid2}')
        # Remove both directions
        self.neighborhood[nid1].remove(nid2)
        self.neighborhood[nid2].remove(nid1)
        if (nid1, nid2) in self.interaction_metadata.keys():
            del self.interaction_metadata[(nid1, nid2)]
        if (nid2, nid1) in self.interaction_metadata.keys():
            del self.interaction_metadata[(nid2, nid1)]

    def remove_node(self, nid):
        if nid not in self.neighborhood.keys():
            raise KeyError(f'Invalid node id: {nid}')
        # remove all interactions with nid
        for other_nid in self.neighborhood[nid]:
            self.neighborhood[other_nid].remove(nid)
            if (nid, other_nid) in self.interaction_metadata.keys():
                del self.interaction_metadata[(nid, other_nid)]
            if (other_nid, nid) in self.interaction_metadata.keys():
                del self.interaction_metadata[(other_nid, nid)]
        if nid in self.node_metadata.keys():
            del self.node_metadata[nid]
        del self.neighborhood[nid]

    def get_neighbors(self, nid) -> set:
        neighbors = self.neighborhood.get(nid)
        if neighbors is not None:
            return neighbors
        return set([])

    def get_node_metadata(self, nid) -> set:
        return self.node_metadata.get(nid)

    def get_interaction_metadata(self, nid1, nid2) -> set:
        return self.interaction_metadata.get((nid1, nid2))

    def __repr__(self):
        return '\n'.join([f'{k}: ' + ' '.join([f'({k}, {v})' for v in v_set]) for k, v_set in self.neighborhood.items()])


def linear_to_pairwise_neighborhood(Np:int, boundary_adjustment=False) -> PairwiseNeighborhood:
    neighborhood = dict()
    # Special case if only 1 node
    if Np == 1:
        neighborhood[0] = set([])            
        return PairwiseNeighborhood(neighborhood)
    neighborhood[0] = set()
    for nid in range(Np-1):
        neighborhood[nid].add(nid+1)
        neighborhood[nid+1] = set([nid])
    if boundary_adjustment:
        neighborhood[0].add(Np-1)
        neighborhood[Np-1].add(0)
    return PairwiseNeighborhood(neighborhood)


def rect2d_to_pairwise_von_neumann_neighborhood(rows, cols, boundary_adjustment=False) -> PairwiseNeighborhood:
    neighborhood = dict()
    # Count from top left -> right
    # For example:
    #  0  1  2  3
    #  4  5  6  7
    #  8  9 10 ..
    # Given (5), neighbors are 1, 4, 5, 9
    neighborhood = dict()
    i = 0
    # first row
    neighborhood[0] = set()
    for _ in range(cols-1):
        neighborhood[i].add(i+1)
        neighborhood[i+1] = set([i])
        i += 1
    # second row onward 
    for m in range(1, rows):
        i += 1
        neighborhood[i] = set()
        for _ in range(cols-1):
            # top
            neighborhood[i].add(i-cols)
            neighborhood[i-cols].add(i)
            # right
            neighborhood[i].add(i+1)
            neighborhood[i+1] = set([i])
            i += 1
        # last on row
        # top
        neighborhood[i].add(i-cols)
        neighborhood[i-cols].add(i)
    
    if boundary_adjustment:
        # connect top and bottom
        top_row_ids = [i for i in range(0, cols)]
        bottom_row_ids = [i for i in range((rows-1)*cols, rows*cols)]
        for nid1, nid2 in zip(top_row_ids, bottom_row_ids):
            neighborhood[nid1].add(nid2)
            neighborhood[nid2].add(nid1)
        # connect sides
        for m in range(rows):
            left_nid = m * cols
            right_nid = (m * cols) + (cols - 1)
            neighborhood[left_nid].add(right_nid)
            neighborhood[right_nid].add(left_nid)

    return PairwiseNeighborhood(neighborhood)


def rect2d_to_pairwise_moore_neighborhood(rows, cols, boundary_adjustment=False) -> PairwiseNeighborhood:
    # Count from top left -> right
    # For example:
    #  0  1  2  3
    #  4  5  6  7
    #  8  9 10 ..
    # Given (5), neighbors are 0, 1, 2, 4, 6, 8, 9, 10

    return NotImplementedError()

    # neighborhood = rect2d_to_pairwise_von_neumann_neighborhood(rows, cols, boundary_adjustment=boundary_adjustment)
    # # add diagonal neighbors
    # for m in range(1, rows):
    #     left_nid = m * cols
    #     top_right_nid = left_nid-(cols-1)
    #     neighborhood.neighborhood[left_nid].add(top_right_nid)
    #     neighborhood.neighborhood[top_right_nid].add(left_nid)

    #     for n in range(1, cols-1):
    #         nid = left_nid + n
    #         top_left_nid = nid-(cols+1)
    #         neighborhood.neighborhood[nid].add(top_left_nid)
    #         neighborhood.neighborhood[top_left_nid].add(nid)

    #         top_right_nid = nid-(cols-1)
    #         neighborhood.neighborhood[nid].add(top_right_nid)
    #         neighborhood.neighborhood[top_right_nid].add(nid)

    #     right_nid = (m * cols) + (cols - 1)
    #     top_left_nid = right_nid - (cols + 1)
    #     neighborhood.neighborhood[right_nid].add(top_left_nid)
    #     neighborhood.neighborhood[top_left_nid].add(right_nid)

    # if boundary_adjustment:
    #     # connect top and bottom diagonals
    #     top_row_ids = [i for i in range(0, cols)]
    #     bottom_row_ids = [i for i in range((rows-1)*cols, rows*cols)]
    #     for i, nid in enumerate(top_row_ids[:-1]):
    #         # top-left
    #         neighborhood.neighborhood[nid].add(bottom_row_ids[i-1])
    #         neighborhood.neighborhood[bottom_row_ids[i-1]].add(nid)
    #         # top-right
    #         neighborhood.neighborhood[nid].add(bottom_row_ids[i-1])
    #         neighborhood.neighborhood[bottom_row_ids[i-1]].add(nid)
    #     i =+ 1
    #     nid = top_row_ids[i]
    #     # bottom-left
    #     neighborhood.neighborhood[nid].add(bottom_row_ids[i-1])
    #     neighborhood.neighborhood[bottom_row_ids[i-1]].add(nid)
    #     # bottom-right
    #     neighborhood.neighborhood[nid].add(bottom_row_ids[0])
    #     neighborhood.neighborhood[bottom_row_ids[0]].add(nid)
    #     # connect side diagonals
    #     for m in range(1, rows):
    #         left_nid = m * cols
    #         right_nid = (m * cols) + (cols - 1)
    #         # top-left
    #         neighborhood.neighborhood[left_nid].add(left_nid-1)
    #         neighborhood.neighborhood[left_nid-1].add(left_nid)
    #         # top-right
    #         neighborhood.neighborhood[right_nid].add(left_nid+1-cols)
    #         neighborhood.neighborhood[left_nid+1-cols].add(right_nid)

    # return neighborhood


def hex2d_to_pairwise_neighborhood(rows, cols, boundary_adjustment=False) -> PairwiseNeighborhood:
    neighborhood = dict()
    # Count from top left -> right
    # Odd cols are 1/2 right offset
    # For example:
    # rows=4, cols=4
    #  0   1   2   3    [4]
    #    4   5   6   7  [4]
    #  8   9  10   11   [4]
    #   12  13  14  15  [4]
    i = 0
    # first row
    neighborhood[0] = set()
    for _ in range(cols-1):
        neighborhood[i].add(i+1)
        neighborhood[i+1] = set([i])
        i += 1
    
    # second row onwards
    for m in range(1, rows):
        i += 1
        neighborhood[i] = set()

        top_left_offset = cols if m % 2 == 1 else cols+1
        # even row 2, 4, 6, ...
        if m % 2 == 0:
            # top-right
            neighborhood[i].add(i-top_left_offset+1)
            neighborhood[i-top_left_offset+1].add(i)
            # right
            neighborhood[i].add(i+1)
            neighborhood[i+1] = set([i])
            i += 1

        start_pos = 0 if m % 2 == 1 else 1
        end_pos = cols-1
        for _ in range(start_pos, end_pos):
            # top-left
            neighborhood[i].add(i-top_left_offset)
            neighborhood[i-top_left_offset].add(i)
            # top-right
            neighborhood[i].add(i-top_left_offset+1)
            neighborhood[i-top_left_offset+1].add(i)
            # right
            neighborhood[i].add(i+1)
            neighborhood[i+1] = set([i])
            i += 1
        # top-left
        neighborhood[i].add(i-top_left_offset)
        neighborhood[i-top_left_offset].add(i)
        # even row 2, 4, 6, ...
        if m % 2 == 0:
            # top-right
            neighborhood[i].add(i-top_left_offset+1)
            neighborhood[i-top_left_offset+1].add(i)

    if boundary_adjustment:
        if rows % 2 != 0:
            raise Exception(f'Cannot perform boundary adjustment with odd number of rows: {rows}')
        top_row_ids = [i for i in range(0, cols)]
        bottom_row_ids = [i for i in range((rows-1)*cols, rows*cols)]
        # Top to bottom
        for i, nid in enumerate(top_row_ids):
            # top-left
            top_left_nid = bottom_row_ids[i-1]
            neighborhood[nid].add(top_left_nid)
            neighborhood[top_left_nid].add(nid)
            # top-right
            top_right_nid = bottom_row_ids[i]
            neighborhood[nid].add(top_right_nid)
            neighborhood[top_right_nid].add(nid)
        # Top row sides
        neighborhood[0].add(cols-1)
        neighborhood[cols-1].add(0)
        # rows
        for r in range(1, rows):
            left_nid = r * cols
            right_nid = (r * cols) + (cols - 1)
            # sides
            neighborhood[left_nid].add(right_nid)
            neighborhood[right_nid].add(left_nid)
            # offsets
            if r % 2 == 1:
                top_right_nid = (r-1) * cols
                neighborhood[right_nid].add(top_right_nid)
                neighborhood[top_right_nid].add(right_nid)
            else:
                top_left_nid = left_nid - 1
                neighborhood[left_nid].add(top_left_nid)
                neighborhood[top_left_nid].add(left_nid)

    return PairwiseNeighborhood(neighborhood)


def rect3d_to_pairwise_neighborhood(Np:int, boundary_adjustment=False) -> PairwiseNeighborhood:
    return NotImplementedError()
    # neighborhood = dict()
    # return PairwiseNeighborhood(neighborhood)


def hex3d_to_pairwise_neighborhood(Np:int, boundary_adjustment=False) -> PairwiseNeighborhood:
    return NotImplementedError()
    # neighborhood = dict()
    # return PairwiseNeighborhood(neighborhood)
