import numpy as np
import pdb

def find_order(refined):
    """
    Find the index array to use to sort the ``refined`` and ``density`` arrays.
    """

    order = np.zeros(refined.shape)

    if not refined[0]:
        return [0]


    def find_nested(i):
        cells = [i]
        for cell in range(8):
            i += 1

            if i >= len(refined):
                return i,np.hstack(cells)

            if refined[i]:
                parent = i
                i, sub_cells = find_nested(i)
                cells.append(sub_cells)
            else:
                cells.append(i)
        cells = [cells[j] for j in [0,1,5,3,7,2,6,4,8]]
        return i, np.hstack(cells)
    
    return find_nested(0)[1]
    
