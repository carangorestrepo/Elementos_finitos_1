function element_num = grid_ql_element_num ( nelemx, nelemy )

%*****************************************************************************80
%
%% grid_ql_element_num() counts the elements in a grid of QL quadrilaterals.
%
%  Example:
%
%    Input:
%
%      NELEMX = 3, NELEMY = 2
%
%    Output:
%
%      ELEMENT_NUM = NELEMX * NELEMY = 6
%
%  Licensing:
%
%    This code is distributed under the MIT license.
%
%  Modified:
%
%    20 June 2005
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    integer NELEMX, NELEMY, the number of elements along the
%    X and Y directions.  The number of elements generated will be
%    NELEMX * NELEMY.
%
%  Output:
%
%    integer ELEMENT_NUM, the number of elements in the grid.
%
  element_num = nelemx * nelemy;

  return
end
