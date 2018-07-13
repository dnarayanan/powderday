from __future__ import print_function

import numpy as np


def hyperion_octree_stats(refined):

    if not (len(refined) - 1) % 8 == 0:
        raise ValueError("refined should have shape 8 * n + 1")

    refined = np.array(refined, dtype=bool, copy=False)

    print("Number of items in refined            : {0}".format(len(refined)))
    print("Number of True values in refined      : {0}".format(np.sum(refined)))
    print("Number of False values in refined     : {0}".format(np.sum(~refined)))

    def check_recursive(refined, current_i=0, max_level=0):
        
  

        if refined[current_i]:
            current_i += 1
            max_levels = []
            for i in range(8):
                current_i, max_level_indiv = check_recursive(refined, current_i, max_level+1)
                max_levels.append(max_level_indiv)
            max_level = max(max_levels)
        else:
            current_i += 1
        return current_i, max_level

    try:
        final_i, max_level = check_recursive(refined)
    except IndexError as e:
        max_level = 'unknown'
        consistent = False
    else:
        consistent = True

    print("Array is self-consistent for Hyperion : {0}".format("yes" if consistent else "no"))
    print("Maximum number of levels              : {0}".format(max_level))

    return max_level
if __name__ == "__main__":
    hyperion_octree_stats([True, False, False, False, False, False, False, False, False])
