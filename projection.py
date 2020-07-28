import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os


def getNVP(key):
    # return normal vector and a point on the plane
    # in this case they are the same
    if key == '1':
        nvp = np.array([0, 0, 1])
    elif key == '2':
        nvp = np.array([1, 0, 0])
    elif key == '3':
        nvp = np.array([0, 1, 0])
    elif key == '4':
        nvp = np.array([-1, 0, 0])
    elif key == '5':
        nvp = np.array([0, -1, 0])
    else:
        raise Exception("No such key in getNVP!")
    return nvp


def checkLegal(key, p):
    # check if the interaction point is legal
    x = p[0]
    y = p[1]
    z = p[2]
    legal = False
    if key == '1':
        if x <= 1 and x >= -1 and y <= 1 and y >= -1 and z == 1:
            legal = True
    else:
        if z <= 1:
            legal = True
    '''
    elif key == '2':
        if x == 1 and y <= 1 and y >= -1 and z >= 0 and z <= 1:
            legal = True
    elif key == '3':
        if x <= 1 and x >= -1 and y == 1 and z >= 0 and z <= 1:
            legal = True
    elif key == '4':
        if x == -1 and y <= 1 and y >= -1 and z >= 0 and z <= 1:
            legal = True
    elif key == '5':
        if x <= 1 and x >= -1 and y == -1 and z >= 0 and z <= 1:
            legal = True
    else:
        raise Exception("No such key in checkLegal!")
    '''
    return legal


def getInter(key, p):
    # get intersection between line from origin to point p and key plane
    # legal => 1; => 0 (op is parallel to key plane); => -1 (op's inter has negative z)
    legal = 1
    inter = np.zeros(3)

    l = p
    n = getNVP(key)
    p0 = n

    denom = l @ n
    if denom == 0:
        legal = 0
    else:
        t = (p0 @ n) / denom
        inter = t * l
        z = inter[2]
        if z < 0:
            legal = -1
    return inter, legal


def getAlterInter(key, illegalP, legalP):
    # move illegal point along the direction of the legal point
    # and compute the intersection
    # alpha cannot be large to avoid erros in z=1 plane.
    alpha = 0.1
    a0 = 0
    a1 = 1
    iter = 10
    num = 0
    direction = legalP - illegalP
    while num < iter:
        num += 1
        newP = illegalP + alpha * direction
        # should ensure new inter is outside the rectangle
        # if inside, a smaller value should be adopted
        inter, legal = getInter(key, newP)
        if legal == 1:
            if checkLegal(key, inter):
                # inter is inside key rectangle
                a1 = alpha
                alpha = (alpha + a0) / 2
            else:
                break
        else:
            # legal == -1/0
            a0 = alpha
            alpha = (alpha + a1) / 2

    if num == iter:
        raise Exception("Max iteration reached!")
    return inter


def getFinalInters(key, test, inters, legals):
    # get corrected inters
    # sm: num of legal inters
    # flags: record which inter is the altered one
    flags = np.zeros(3)
    sm = 0
    for i in range(3):
        if legals[i] == 1:
            sm += 1
    # direction corresponding to each inter
    # can be optimized (eg. for 1 legal point, only 2 points and dirs are needed)
    # use zMax to filter some obvious illegal projections
    if sm == 0:
        # need to change!! (for cases with 3 z<0 points) + flag
        print("No actual projection detected!")
        # ensure newinters have no intersections with rectangles!
        newinters = np.ones((3, 3)) * (-10)
    elif sm == 1:
        newinters = np.zeros((3, 3))
        for i in range(3):
            if legals[i] < 1:
                for j in range(2):
                    k = (i + j + 1) % 3
                    if legals[k] == 1:
                        legalP = test[k]
                        illegalP = test[i]
                        altInter = getAlterInter(key, illegalP, legalP)
                        newinters[i] = altInter
                        flags[i] = 1
            else:
                newinters[i] = inters[i]
    elif sm == 2:
        newinters = np.zeros((4, 3))
        flags = np.zeros(4)
        for i in range(3):
            j = 0
            for i in range(3):
                if legals[i] < 1:
                    illegalP = test[i]
                    legalP = test[(i + 2) % 3]
                    altInter = getAlterInter(key, illegalP, legalP)
                    newinters[j] = altInter
                    flags[j] = 1
                    j += 1

                    legalP = test[(i + 1) % 3]
                    altInter = getAlterInter(key, illegalP, legalP)
                    newinters[j] = altInter
                    flags[j] = 1
                    j += 1
                else:
                    newinters[j] = inters[i]
                    j += 1
    elif sm == 3:
        newinters = inters
    else:
        raise Exception("Sum of legals out of bound!")
    return newinters, flags


def getInters(key, test):
    # get insections between triangle test and key plane
    # inters is in shape (3,3) or (4,3)
    inters = np.zeros((3, 3))
    legals = np.zeros(3)
    for i in range(3):
        inter, legal = getInter(key, test[i])
        inters[i] = inter
        legals[i] = legal

    inters, flags = getFinalInters(key, test, inters, legals)
    return inters, flags


def getInterRec(rec1, rec2):
    # get the interaction region between two rectangles
    # each rectangle is in form rec1 = array[[xmin, ymin], [xmax, ymax]]
    exist = True
    recInter = np.zeros((2, 2))
    startX = max(rec1[0][0], rec2[0][0])
    endX = min(rec1[1][0], rec2[1][0])
    startY = max(rec1[0][1], rec2[0][1])
    endY = min(rec1[1][1], rec2[1][1])
    if startX >= endX or startY >= endY:
        exist = False
    else:
        recInter[0, 0] = startX
        recInter[0, 1] = startY
        recInter[1, 0] = endX
        recInter[1, 1] = endY
    return recInter, exist


def getInterRegion(key, inters):
    # get intersection region between inters and the key plane
    # return recInter is in shape (2,2) => 2D
    p1 = np.min(inters, axis=0).reshape(1, -1)
    p2 = np.max(inters, axis=0).reshape(1, -1)
    rec = np.concatenate((p1, p2), axis=0)
    rec2 = g2l(key, rec)

    if key == '1':
        rec1 = np.array([[-1, -1], [1, 1]])
    elif key == '2':
        rec1 = np.array([[-1, 0], [1, 1]])
    elif key == '3':
        rec1 = np.array([[-1, 0], [1, 1]])
    elif key == '4':
        rec1 = np.array([[-1, 0], [1, 1]])
    elif key == '5':
        rec1 = np.array([[-1, 0], [1, 1]])
    else:
        raise Exception("No such key in getInterRegion")
    recInter, exist = getInterRec(rec1, rec2)
    return recInter, exist


def g2l(key, p):
    # from global coordinates to local
    # 3D => 2D
    # p is in shape (-1,3); pl is in shape (-1,2)
    p = p.reshape(-1, 3)
    if key == '1':
        pl = p[:, (0, 1)]
    elif key == '2':
        pl = p[:, (1, 2)]
    elif key == '3':
        pl = p[:, (0, 2)]
    elif key == '4':
        pl = p[:, (1, 2)]
    elif key == '5':
        pl = p[:, (0, 2)]
    else:
        raise Exception("No such key in g2l!")
    return pl


def l2g(key, p):
    # from local coordinates to global
    # 2D => 3D
    # p is in shape (-1,2); pg is in shape (-1,3)
    p = p.reshape(-1, 2)
    m = p.shape[0]

    if key == '1':
        pg = np.concatenate((p[:, 0].reshape(-1, 1), p[:, 1].reshape(-1, 1),
                             (1 * np.ones(m)).reshape(-1, 1)),
                            axis=1)
    elif key == '2':
        pg = np.concatenate(((1 * np.ones(m)).reshape(-1, 1), p[:, 0].reshape(
            -1, 1), p[:, 1].reshape(-1, 1)),
                            axis=1)
    elif key == '3':
        pg = np.concatenate(
            (p[:, 0].reshape(-1, 1),
             (1 * np.ones(m)).reshape(-1, 1), p[:, 1].reshape(-1, 1)),
            axis=1)
    elif key == '4':
        pg = np.concatenate(((-1 * np.ones(m)).reshape(-1, 1), p[:, 0].reshape(
            -1, 1), p[:, 1].reshape(-1, 1)),
                            axis=1)
    elif key == '5':
        pg = np.concatenate(
            (p[:, 0].reshape(-1, 1),
             (-1 * np.ones(m)).reshape(-1, 1), p[:, 1].reshape(-1, 1)),
            axis=1)
    else:
        raise Exception("No such key in l2g!")
    return pg


def coo2mesh(key, p):
    # from local coordinates to mesh indexes
    # p is a single point in shape (,2)
    index = np.zeros(2)
    if key == '1':
        index[1] = min(int((p[0] + 1) * n), 2 * n - 1)
        index[0] = min(int((p[1] + 1) * n), 2 * n - 1)
    else:
        index[1] = min(int((p[0] + 1) * n), 2 * n - 1)
        index[0] = min(int(p[1] * n), n - 1)
    return index


def mesh2coo(key, index):
    # from mesh indexes to local coordinates of the center of this grid
    # index is indexes of a point in shape (,2)
    p = np.zeros(2)
    if key == '1':
        p[0] = -1 + 1 / (2 * n) + index[1] / n
        p[1] = -1 + 1 / (2 * n) + index[0] / n
    else:
        p[0] = -1 + 1 / (2 * n) + index[1] / n
        p[1] = 1 / (2 * n) + index[0] / n
    return p


def getDir(inters):
    # get the direction of each edge
    # inters is in shape (3.3) or (4,3) which forms a convex polygon
    m, n = inters.shape
    dirs = np.zeros((m, n))
    for i in range(m):
        dirs[i] = inters[(i + 1) % m] - inters[i]
    return dirs


def checkIn(p, inters, dirs):
    # check if point p inside the polygon
    isIn = False
    m = inters.shape[0]
    values = []
    for i in range(m):
        d1 = p - inters[i]
        d2 = dirs[i]
        # c has 2 zeros and a value
        c = np.cross(d1, d2)
        values.append(np.sum(c))
    if min(values) >= 0 or max(values) <= 0:
        isIn = True
    return isIn


n = 10
mesh = {
    '1': np.zeros((2 * n, 2 * n)),
    '2': np.zeros((n, 2 * n)),
    '3': np.zeros((n, 2 * n)),
    '4': np.zeros((n, 2 * n)),
    '5': np.zeros((n, 2 * n))
}

#test = np.array([[0,1,2],[4,0,2],[4,4,2]])
#test = np.array([[0,1,2],[0,0,2],[4,4,2]])
#test = np.array([[0,1,2],[0,0,2],[10,2,2]])
#test = np.array([[10,0,2],[-10,-10,2],[-10,10,2]])/15
#test = np.array([[-2,1,2],[4,0,2],[4,4,2]])
#test = np.array([[0, 1, 2], [4, 0, 2], [4, 2, 2]])
test = np.array([[-3, 1, 2], [4, 0, 2], [4, 2, 0]], dtype=np.float64)

# z = 0 should be avoided in this method!
for i in range(3):
    if test[i, 2] == 0:
        test[i, 2] = 1e-6

for key in mesh.keys():
    print(key)
    inters, flags = getInters(key, test)

    dirs = getDir(inters)
    #print(inters)
    #print(dirs)

    recInter, exist = getInterRegion(key, inters)
    #print(recInter, exist)
    if exist:
        minIndex = coo2mesh(key, recInter[0])
        maxIndex = coo2mesh(key, recInter[1])
        #print(minIndex, maxIndex)
        for i in np.arange(minIndex[0], maxIndex[0] + 1):
            for j in np.arange(minIndex[1], maxIndex[1] + 1):
                index = np.array([i, j])
                p = mesh2coo(key, index)
                p = l2g(key, p)
                if checkIn(p, inters, dirs):
                    #print(i,j)
                    mesh[key][int(i), int(j)] = 1
                    #print(p)
        print(np.sum(mesh[key]))


def draw(x, y, z, dx, dy, dz, tri, mesh, color='blue'):
    fig = plt.figure(figsize=(8, 8))
    ax = Axes3D(fig)
    xx = [x, x, x + dx, x + dx, x]
    yy = [y, y + dy, y + dy, y, y]
    kwargs = {'alpha': 0.5, 'color': color}
    # draw lines to form face
    ax.plot3D(xx, yy, [z] * 5, **kwargs)
    ax.plot3D(xx, yy, [z + dz] * 5, **kwargs)
    ax.plot3D([x, x], [y, y], [z, z + dz], **kwargs)
    ax.plot3D([x, x], [y + dy, y + dy], [z, z + dz], **kwargs)
    ax.plot3D([x + dx, x + dx], [y + dy, y + dy], [z, z + dz], **kwargs)
    ax.plot3D([x + dx, x + dx], [y, y], [z, z + dz], **kwargs)

    tri = tri.T
    Xt = tri[0]
    Yt = tri[1]
    Zt = tri[2]
    ax.plot_trisurf(Xt, Yt, Zt, **kwargs)

    for i in range(3):
        X = [0, tri[0, i]]
        Y = [0, tri[1, i]]
        Z = [0, tri[2, i]]
        ax.plot3D(X, Y, Z, **kwargs)

    for key in mesh.keys():
        meshes = mesh[key]
        m, n = meshes.shape
        for i in range(m):
            for j in range(n):
                index = np.array([i, j])
                p = mesh2coo(key, index)
                p = l2g(key, p)
                if meshes[i, j]:
                    ax.scatter(p[0, 0], p[0, 1], p[0, 2], c='r')
                else:
                    pass
                    #ax.scatter(p[0,0],p[0,1],p[0,2],c='k')
    #ax.set_aspect('equal')
    xmin = min(-1, min(Xt))
    xmax = max(1, max(Xt))
    ymin = min(-1, min(Yt))
    ymax = max(1, max(Yt))
    zmin = min(0, min(Zt))
    zmax = max(1, max(Zt))
    xl = xmax - xmin
    yl = ymax - ymin
    zl = zmax - zmin
    l = max(xl, yl, zl)
    ax.set_xlim(xmin - 0.5, xmin + l + 0.5)
    ax.set_ylim(ymin - 0.5, ymin + l + 0.5)
    ax.set_zlim(zmin - 0.5, zmin + l + 0.5)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()


draw(-1, -1, 0, 2, 2, 1, test, mesh)
