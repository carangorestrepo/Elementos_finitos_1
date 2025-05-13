function [ element_neighbor ] = grid_t3_neighbor ( nelemx, nelemy )

%*****************************************************************************80
%
%% grid_t3_neighbor() produces the element neighbor array.
%
%  Discussion:
%
%    A rectangular grid is assumed, which can be thought of as a
%    grid of (nelemx+1) x (nelemy+1) evenly spaced points, or as
%    a grid of nelemx x nelemy squares, each square being subdivided
%    into a pair of triangles.
%
%    Node numbering starts at the lower left, moves right first, and
%    then upwards.
%
%    Element numbering goes similarly, square by square, with the
%    lower triangle of each pair numbered first.
%
%    The element node array records, for each element, the indices of
%    the three nodes that form the vertices of the triangle, in counter
%    clockwise order, starting with the node associated with the 90 degree
%    angle.
%
%    The element neighbor array records the index of the element opposite
%    each node.
%
%  Example:
%
%    NELEMX = 3, NELEMY = 2
%
%    element_node =
%       1,  2,  5;
%       6,  5,  2;
%       2,  3,  6;
%       7,  6,  3;
%       3,  4,  7;
%       8,  7,  4;
%       5,  6,  9;
%      10,  9,  6;
%       6,  7, 10;
%      11, 10,  7;
%       7,  8, 11;
%      12, 11,  8.
%
%    element_neighbor =
%       2, -1, -1;
%       1,  3,  7;
%       4,  2, -1;
%       3,  5,  9;
%       6,  4, -1;
%       5, -1, 11;
%       8, -1,  2;
%       7,  9, -1;
%      10,  8,  4;
%       9, 11, -1;
%      12, 10,  6;
%      11, -1, -1;
%
%  Grid:
%
%    9---10---11---12     +----+----+----+     +----+----+----+
%    |\   |\   |\   |     |\ 8 |\10 |\12 |     |\ * |\ * |\ * |
%    | \  | \  | \  |     | \  | \  | \  |     | \79| \9E| \E*|
%    |  \ |  \ |  \ |     |  \ |  \ |  \ |     |*8\ |8X\ |XT\ |
%    |   \|   \|   \|     |  7\|  9\| 11\|     | 2 \| 4 \| 6 \|
%    5----6----7----8     +----+----+----+     +----+----+----+
%    |\   |\   |\   |     |\ 2 |\ 4 |\ 6 |     |\ 7 |\ 9 |\ E |
%    | \  | \  | \  |     | \  | \  | \  |     | \13| \35| \5*|
%    |  \ |  \ |  \ |     |  \ |  \ |  \ |     |*2\ |24\ |46\ |
%    |   \|   \|   \|     |  1\|  3\|  5\|     | * \| * \| * \|
%    1----2----3----4     +----+----+----+     +----+----+----+
%      Node indexing      Element indexing     Neighbor indexing
%
%  Reference Element T3:
%
%    |
%    1  3
%    |  |\
%    |  | \
%    S  |  \
%    |  |   \
%    |  |    \
%    0  1-----2
%    |
%    +--0--R--1-->
%
%  Licensing:
%
%    This code is distributed under the MIT license.
%
%  Modified:
%
%    06 April 2019
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    integer NELEMX, NELEMY, the number of elements along the
%    X and Y directions.  The number of elements generated will be
%    2 * NELEMX * NELEMY.
%
%  Output:
%
%    integer ELEMENT_NEIGHBOR(3,2*NELEMX*NELEMY), the indices
%    of neighbor elements.
%

  e = 0;

  for j = 1 : nelemy
    for i = 1 : nelemx

      e = e + 1;

      e1 = e + 1;
      if ( i == 1 )
        e2 = -1;
      else
        e2 = e - 1;
      end
      if ( j == 1 )
        e3 = -1;
      else
        e3 = e - 2 * nelemx + 1;
      end

      element_neighbor(1,e) = e1;
      element_neighbor(2,e) = e2;
      element_neighbor(3,e) = e3;

      e = e + 1;

      e1 = e - 1;
      if ( i == nelemx )
        e2 = -1;
      else
        e2 = e + 1;
      end
      if ( j == nelemy )
        e3 = -1;
      else
        e3 = e + 2 * nelemx - 1;
      end

      element_neighbor(1,e) = e1;
      element_neighbor(2,e) = e2;
      element_neighbor(3,e) = e3;

    end
  end

  return
end
