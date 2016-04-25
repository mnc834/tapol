%%%-------------------------------------------------------------------
%%% @author PKQ874
%%% @copyright (C) 2016, <COMPANY>
%%% @doc
%%%
%%% @end
%%% Created : 14. мар 2016 18:31
%%% The module is a set of functions that calculate least square solution
%%% for a given array of pairs {X, Y}. The solution is a polinomial of an
%%% abitary polinomial degree, presented in a form of list of floats.
%%% The list of floats is the least of the polinomial coefficiants when the
%%% coefficient with highest order comes first
%%%-------------------------------------------------------------------
-module(tapol_lss).
-author("PKQ874").

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

-include("tapol.hrl").

-type e_line() :: {[float()], float()}.
%% @doc equuasion line, left-hand side coefficients and right-hand side member

-record(matrix, {
  lines :: [e_line()],                  %%lines of the equasion matrix
  column_order :: [pos_integer()]       %%list of one-based column numbers
}).


%% API
-export([get_least_squares_solution/2]).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-spec extract_submatrix(#matrix{}) -> #matrix{}.
%% @doc removes the first line and the first column from the matrix
extract_submatrix(M) ->
  remove_first_column(remove_first_line(M)).

-spec find_column_with_lowest_rank(#matrix{}) -> Result
  when Result :: {ok, Colunm_id :: pos_integer(), Column_pos :: pos_integer()} |
                 {error, not_found}.
%% @doc looks for the column with the bigest number of elements that are less than ?EPSILON
%%      returns {ok, column_ID, Zero_based_column_position} if found
%%      or {error, not_found} if there are no elements with less than ?EPSILON in the table
find_column_with_lowest_rank(#matrix{lines = Lines, column_order = Order}) ->
  N = length(Order),
  Line_fun =
    fun({V, C}) ->
      case abs(V) < ?EPSILON of
        true ->
          C + 1;
        false ->
          C
      end
    end,
  Matrix_fun =
    fun({Lhs, _Rhs}, Acc) ->
      Zip_line = lists:zip(Lhs, Acc),
      lists:map(Line_fun, Zip_line)
    end,
  Count_list =
    lists:foldl(Matrix_fun, lists:map(fun(_) -> 0 end, lists:seq(1, N)), Lines),
  Column_list = lists:zip(Order, lists:seq(0, N - 1)),
  Processed_list = lists:zip(Count_list, Column_list),
  Sort_fun =
    fun({V1, _}, {V2, _}) ->
      (V1 > V2)
    end,
  [{Sum, {Column_id, Column_pos}} | _] = lists:sort(Sort_fun, Processed_list),
  case Sum =:= 0 of
    true ->
      {error, not_found};
    false ->
      {ok, Column_id, Column_pos}
  end.

-spec move_column_to_first_place(#matrix{}, Column_id :: pos_integer()) -> {ok, #matrix{}} | {error, Reason :: term()}.
%% @doc makes a column with ID Column_id to be the first column
%%      returns {ok, #matrix{}} in case of success of {error, Reason} in case of failure
move_column_to_first_place(#matrix{lines = Lines, column_order = Columns} = M, Column_id) ->
  Pos_fun =
    fun(Column, {not_found, I}) when Column =:= Column_id ->
      {found, I};
      (_, {not_found, I}) ->
        {not_found, I + 1};
      (_, {found, _} = Res) ->
        Res
    end,
  case lists:foldl(Pos_fun, {not_found, 0}, Columns) of
    {found, 0} ->
      %%nothing to be changes since the first element should remain the first
      {ok, M};
    {found, Pos} ->
      %%for every Line including Columns move the Pos + 1 element to become first
      Move_fun =
        fun Fun({Lhs, Rhs}) ->
          {Fun(Lhs), Rhs};
          Fun(L) when is_list(L)->
            {L1, [Head | Tail]} = lists:split(Pos, L),
            [Head | L1] ++ Tail
        end,
      New_lines = lists:map(Move_fun, Lines),
      New_columns = Move_fun(Columns),
      {ok, M#matrix{lines = New_lines, column_order = New_columns}};
    {not_found, _} ->
      {error, {column_id_not_present, Column_id}}
  end.

-spec sort_lines(#matrix{}) -> #matrix{}.
%% @doc sorts the matrix lines based on the absolute value of the first element
%%      puts the line with the biggest absolute first element to the first position
sort_lines(#matrix{lines = Lines} = M) ->
  Sort_fun =
    fun({[A | _], _Rhs_a}, {[B | _], _Rhs_b}) ->
      abs(A) > abs(B)
    end,
  New_lines = lists:sort(Sort_fun, Lines),
  M#matrix{lines = New_lines}.

-spec gauss_elimination(#matrix{}) -> {ok, #matrix{}, pos_integer()} | {error, Reason :: term()}.
%% @doc find the column with lowest rank, moves it to the first place, sorts the matrix,
%%      eliminates all members in the first column except for the very first one and
%%      returns the updated matrix and the initial zero-based position on the column
%%      with lowest rank
%%      in case when all elements in the column with lowest rank are less than ?EPSILONE
%%      returns {error, no_simple_solution}
gauss_elimination(#matrix{} = M) ->
  {Permut_m, Pos} =
    case find_column_with_lowest_rank(M) of
      {ok, Id, Pos_val} ->
        {ok, M1} = move_column_to_first_place(M, Id),
        {M1, Pos_val};
      {error, not_found} ->
        %%No columns to be permute
        {M, 0}
    end,
  Ready_m = sort_lines(Permut_m),
  Normalize_line =
    fun({[Lhs_h | _] = Lhs, Rhs}) ->
      {[V / Lhs_h || V <- Lhs], Rhs / Lhs_h}
    end,
  Subtract_lines =
    fun({Lhs_a, Rhs_a}, {Lhs_b, Rhs_b}) ->
      Lhs =
        lists:map(fun({A, B}) -> A - B end, lists:zip(Lhs_a, Lhs_b)),
      Rhs = Rhs_a - Rhs_b,
      {Lhs, Rhs}
    end,
  #matrix{lines = [{[Top_left_element | _], _} = Line_1 | Rest_lines]} = Ready_m,
  %%checking if the elimination is possible
  case abs(Top_left_element) > ?EPSILON of
    true ->
      Norm_line_1 = Normalize_line(Line_1),
      %%eliminating the first elements in all lines
      Eliminate_fun =
        fun({[E1 | _], _} = Line, true)
          when abs(E1) > ?EPSILON ->
          Norm_line = Normalize_line(Line),
          {Subtract_lines(Norm_line, Norm_line_1), true};
          (Line, true) ->
            %%the first element of the line is less than ?EPSILONE
            %%stop processing
            {Line, false};
          (Line, false) ->
            %%processing is stopped
            {Line, false}
        end,
      {New_rest_lines, _Flag} =
        lists:mapfoldl(Eliminate_fun, true, Rest_lines),
      New_matrix = Ready_m#matrix{lines = [Norm_line_1 | New_rest_lines]},
      {ok, New_matrix, Pos};
    false ->
      {error, no_simple_solution}
  end.

-spec gauss_solution(#matrix{}) -> {ok, [float()]} | {error, Reason}
  when Reason :: no_simple_solution.
%% @doc solves linear equqsions system, returns {ok, solution vector} of
%%      {error, Reason} in case of solution cannot be found
gauss_solution(#matrix{lines = [{[V], Rhs}]}) ->
  %%base case
  case abs(V) > ?EPSILON of
    true ->
      {ok, [Rhs / V]};
    false ->
      {error, no_simple_solution}
  end;
gauss_solution(M) ->
  case gauss_elimination(M) of
    {ok, M1, Pos} ->
      #matrix{lines = [{Lhs, Rhs} | _]} = M1,
      Sub_matrix = extract_submatrix(M1),
      case gauss_solution(Sub_matrix) of
        {ok, Vector} ->
          Vector_1 = [0 | Vector],
          V = Rhs - lists:foldl(fun({V, A}, S) -> S + V * A end, 0, lists:zip(Vector_1, Lhs)),
          %% have to put V to the right place
          {L1, L2} = lists:split(Pos, Vector),
          Solution = L1 ++ [V | L2],
          {ok, Solution};
        {error, Reason} ->
          {error, Reason}
      end;
    {error, Error} ->
      {error, Error}
  end.

-spec get_least_squares_solution([{X, Y}], Polynomial_degree) -> Result
  when X :: float(),
       Y :: float(),
       Polynomial_degree :: pos_integer(),
       Result :: {ok, tapol_epol:e_polynomial()} | {error, Reason :: term()}.
%% @doc calculates polynomial approximation for the set of data [{X, Y}] using least squares method
%%      the resulting polynom should have degree equals Polynomial_degree
get_least_squares_solution(V, P_degree) ->
  X_l = 2 * P_degree + 1,
  Y_l = P_degree + 1,
  Zero_list =
    fun(N) ->
      [0 || _ <- lists:seq(1, N)]
    end,
  Pow_fun =
    fun(X, Acc) ->
      {Acc, Acc * X}
    end,
  Gen_pow_line =
    fun(X) ->
      {L, _} = lists:mapfoldl(Pow_fun, 1, [X || _ <- lists:seq(1, X_l)]),
      L
    end,
  Foldl_fun =
    fun({X, Y}, {C_x, C_y}) ->
      Add_fun =
        fun({A, B}) ->
          A + B
        end,
      Pow_line = Gen_pow_line(X),
      New_c_x = lists:map(Add_fun, lists:zip(Pow_line, C_x)),
      {Y_pow_line, _} = lists:split(P_degree + 1, Pow_line),
      Y_line = [Y * Y_pow_line_v || Y_pow_line_v <- Y_pow_line],
      New_c_y = lists:map(Add_fun, lists:zip(Y_line, C_y)),
      {New_c_x, New_c_y}
    end,
  {C_x, C_y} =
    lists:foldl(Foldl_fun, {Zero_list(X_l), Zero_list(Y_l)}, V),

  %%creating the Lhs of the lines
  Get_lhs_fun =
    fun(_X, [_H | T] = L) ->
      {Lhs, _} = lists:split(P_degree + 1, L),
      {Lhs, T}
    end,
  {Lhs_list, _} =
    lists:mapfoldl(Get_lhs_fun, C_x, lists:seq(1, P_degree +1)),
  %%creating matrix lines
  Lines = lists:zip(Lhs_list, C_y),
  %%creating the matrix
  M = #matrix{lines = Lines, column_order = lists:seq(1, P_degree + 1)},

  %Getting solution
  case gauss_solution(M) of
    {ok, R_solution} ->
      %%higher order coefficients should come first
      {ok, lists:reverse(R_solution)};
    Error ->
      Error
  end.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Local funs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-spec remove_first_line(#matrix{}) -> #matrix{}.
%% @doc removes the first line from the matrix
remove_first_line(#matrix{lines = [_H_line| T_line]} = M) ->
  M#matrix{lines = T_line}.

-spec remove_first_column(#matrix{}) -> #matrix{}.
%% @doc removes the first column from the matrix
remove_first_column(#matrix{lines = Lines, column_order = [_H_order | T_order]} = M) ->
  New_lines = [{T_line, Rhs} || {[_H_line | T_line], Rhs} <- Lines],
  M#matrix{lines = New_lines, column_order = T_order}.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Unit tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-ifdef(TEST).
extract_submatrix_test_() ->

  M = #matrix{
    lines =
      [
        {[1, 2, 3, 4], 10},
        {[5, 6, 7, 8], 20},
        {[9, 10, 11, 2], 30},
        {[13, 14, 15, 16], 40}
      ],
    column_order = [1, 2, 3, 4]
  },

  M1 = #matrix{
    lines =
    [
      {[6, 7, 8], 20},
      {[10, 11, 2], 30},
      {[14, 15, 16], 40}
    ],
    column_order = [2, 3, 4]
  },
  M2 = #matrix{
    lines =
    [
      {[11, 2], 30},
      {[15, 16], 40}
    ],
    column_order = [3, 4]
  },
  M3 = #matrix{
    lines =
    [
      {[16], 40}
    ],
    column_order = [4]
  },
  M4 = #matrix{
    lines =
    [
    ],
    column_order = []
  },

  [
    ?_assertEqual(extract_submatrix(M), M1),
    ?_assertEqual(extract_submatrix(M1), M2),
    ?_assertEqual(extract_submatrix(M2), M3),
    ?_assertEqual(extract_submatrix(M3), M4)
  ].

find_column_with_lowest_rank_test_() ->
  E_val = ?EPSILON / 10,
  M = #matrix{
    lines =
    [
      {[1, 2, 3, 4], 10},
      {[5, 6, 7, 8], 20},
      {[9, 10, 11, 2], 30},
      {[13, 14, 15, 16], 40}
    ],
    column_order = [1, 2, 3, 4]
  },
  M1 = #matrix{
    lines =
    [
      {[-1, 2, 3, 4], 10},
      {[-5, 6, 7, 8], 20},
      {[-9, E_val , 11, 2], 30},
      {[-13, 14, 15, 16], 40}
    ],
    column_order = [1, 5, 3, 4]
  },
  M2 = #matrix{
    lines =
    [
      {[-1, 2, 3, 4], 10},
      {[-5, 6, E_val, 8], 20},
      {[-9, E_val , E_val, 2], 30},
      {[-13, 14, 15, 16], 40}
    ],
    column_order = [1, 5, 6, 4]
  },
  M3 = #matrix{
    lines =
    [
      {[-1, 2, 3, 4], 10},
      {[-5, 6, E_val, E_val], 20},
      {[-9, E_val , E_val, E_val], 30},
      {[-13, 14, 15, E_val], 40}
    ],
    column_order = [4, 5, 6, 1]
  },
  M4 = #matrix{
    lines =
    [
      {[E_val, 2, 3, 4], 10},
      {[-5, 6, E_val, E_val], 20},
      {[-9, E_val , E_val, E_val], 30},
      {[-13, 14, 15, E_val], 40}
    ],
    column_order = [4, 5, 6, 1]
  },
  M5 = #matrix{
    lines =
    [
      {[E_val, 2, 3, 4], 10},
      {[E_val, 6, E_val, E_val], 20},
      {[E_val, E_val , E_val, E_val], 30},
      {[E_val, 14, 15, E_val], 40}
    ],
    column_order = [4, 5, 6, 1]
  },

  [
    ?_assertEqual(find_column_with_lowest_rank(M), {error, not_found}),
    ?_assertEqual(find_column_with_lowest_rank(M1), {ok, 5, 1}),
    ?_assertEqual(find_column_with_lowest_rank(M2), {ok, 6, 2}),
    ?_assertEqual(find_column_with_lowest_rank(M3), {ok, 1, 3}),
    ?_assertEqual(find_column_with_lowest_rank(M4), {ok, 1, 3}),
    ?_assertEqual(find_column_with_lowest_rank(M5), {ok, 4, 0})
  ].

move_column_to_first_place_test_() ->
  M = #matrix{
    lines =
    [
      {[1, 2, 3, 4], 10},
      {[5, 6, 7, 8], 20},
      {[9, 10, 11, 12], 30},
      {[13, 14, 15, 16], 40}
    ],
    column_order = [1, 2, 3, 4]
  },
  M1 = #matrix{
    lines =
    [
      {[2, 1, 3, 4], 10},
      {[6, 5, 7, 8], 20},
      {[10, 9, 11, 12], 30},
      {[14, 13, 15, 16], 40}
    ],
    column_order = [2, 1, 3, 4]
  },
  M2 = #matrix{
    lines =
    [
      {[4, 2, 1, 3], 10},
      {[8, 6, 5, 7], 20},
      {[12, 10, 9, 11], 30},
      {[16, 14, 13, 15], 40}
    ],
    column_order = [4, 2, 1, 3]
  },
  [
    ?_assertEqual(move_column_to_first_place(M, 10), {error, {column_id_not_present, 10}}),
    ?_assertEqual(move_column_to_first_place(M, 1), {ok, M}),
    ?_assertEqual(move_column_to_first_place(M, 2), {ok, M1}),
    ?_assertEqual(move_column_to_first_place(M1, 4), {ok, M2})
  ].

sort_lines_test_() ->
  A1 = #matrix{
    lines =
    [
      {[1, 2, 3, 4], 10},
      {[5, 6, 7, 8], 20},
      {[9, 10, 11, 12], 30},
      {[13, 14, 15, 16], 40}
    ],
    column_order = [1, 2, 3, 4]
  },

  A2 = #matrix{
    lines =
    [
      {[13, 14, 15, 16], 40},
      {[9, 10, 11, 12], 30},
      {[5, 6, 7, 8], 20},
      {[1, 2, 3, 4], 10}
    ],
    column_order = [1, 2, 3, 4]
  },

  B1 = #matrix{
    lines =
    [
      {[1, 2, 3, 4], 10},
      {[-50, 6, 7, 8], 20},
      {[9, 10, 11, 12], 30},
      {[-13, 14, 15, 16], 40}
    ],
    column_order = [1, 2, 3, 4]
  },
  B2 = #matrix{
    lines =
    [
      {[-50, 6, 7, 8], 20},
      {[-13, 14, 15, 16], 40},
      {[9, 10, 11, 12], 30},
      {[1, 2, 3, 4], 10}
    ],
    column_order = [1, 2, 3, 4]
  },
  [
    ?_assertEqual(sort_lines(A1), A2),
    ?_assertEqual(sort_lines(B1), B2)
  ].

gauss_elimination_test_() ->
  Compare_matrix =
    fun(#matrix{lines = Lines_a, column_order = Order_a},
        #matrix{lines = Lines_b, column_order = Order_b}) ->
      case Order_a =:= Order_b of
        true ->
          Compare_lines_fun =
            fun({{Lhs_a, Rhs_a}, {Lhs_b, Rhs_b}}) ->
              case tapol_utils:are_equal(Rhs_a, Rhs_b) of
                true ->
                  lists:all(fun({A, B}) -> tapol_utils:are_equal(A, B) end, lists:zip(Lhs_a, Lhs_b));
                false ->
                  {error, right_hand_sides_different, Rhs_a, Rhs_b}
              end
            end,
            lists:all(Compare_lines_fun, lists:zip(Lines_a, Lines_b));
        false ->
          {error, column_orders_different, Order_a, Order_b}
      end
    end,

  A0 = #matrix{
    lines =
    [
      {[1, 0, 3, 4], 10},
      {[5, 0, 7, 8], 20},
      {[9, 0, 11, 12], 30},
      {[13, 0, 15, 16], 40}
    ],
    column_order = [1, 2, 3, 4]
  },
  A1 = #matrix{
    lines =
    [
      {[1, 2, 3, 4], 10},
      {[-5, 6, 7, 8], 20},
      {[9, 10, 11, 12], 30},
      {[13, 14, 15, 16], 40}
    ],
    column_order = [1, 2, 3, 4]
  },
  {ok, A11, 0} = gauss_elimination(A1),
  B11 = #matrix{
    lines =
    [
      {[13 / 13, 14 / 13, 15 / 13, 16 / 13], 40 / 13},
      {[0, 10 / 9 - 14 / 13, 11 / 9 - 15 / 13, 12 / 9 - 16 / 13], 30 / 9 - 40 / 13},
      {[0, 6 / -5 - 14 / 13, 7 / -5 - 15 / 13, 8 / -5 - 16 / 13], 20 / -5 - 40 / 13},
      {[0, 2 - 14 / 13, 3 - 15 / 13, 4 - 16 / 13], 10 - 40 / 13}
    ],
    column_order = [1, 2, 3, 4]
  },

  A2 = #matrix{
    lines =
    [
      {[1, 2, 3, 4], 10},
      {[5, 0, 7, 8], 20},
      {[9, 0, 11, 12], 30},
      {[13, 14, 15, 16], 40}
    ],
    column_order = [1, 2, 3, 4]
  },
  {ok, A21, 1} = gauss_elimination(A2),
  B21 = #matrix{
    lines =
    [
      {[1, 13 / 14, 15 / 14, 16 / 14], 40 / 14},
      {[0, 1 / 2 - 13 / 14, 3 / 2 - 15 / 14, 4 / 2 - 16 / 14], 10 / 2 - 40 / 14},
      {[0, 9, 11, 12], 30},
      {[0, 5, 7, 8], 20}
    ],
    column_order = [2, 1, 3, 4]
  },

  [
    ?_assertEqual({error, no_simple_solution}, gauss_elimination(A0)),
    ?_assertEqual(true, Compare_matrix(A11, B11)),
    ?_assertEqual(true, Compare_matrix(A21, B21))
  ].

gauss_solution_test_() ->
  Check_solution =
    fun(#matrix{lines = Lines}, V) ->
      Sum_fun =
        fun({A, B}, S) ->
          S + A * B
        end,
      Line_fun =
        fun({Lhs, Rhs}) ->
          Sum = lists:foldl(Sum_fun, 0, lists:zip(Lhs, V)),
          tapol_utils:are_equal(Sum, Rhs)
        end,
      lists:all(Line_fun, Lines)
    end,

  A0 = #matrix{
    lines =
    [
      {[1, 0, 3, 4], 10},
      {[5, 0, 7, 8], 20},
      {[9, 0, 11, 12], 30},
      {[13, 0, 15, 16], 40}
    ],
    column_order = [1, 2, 3, 4]
  },

  A1 = #matrix{
    lines =
    [
      {[1, 2, 3], 1},
      {[5, 6, 7], 1},
      {[3, 2, 0], 1}
    ],
    column_order = [1, 2, 3]
  },

  {ok, V1} = gauss_solution(A1),

  A2 = #matrix{
    lines =
    [
      {[1, 2, 2, 2], 10},
      {[4, 5, 6, 7], 20},
      {[1, 4, 1, 4], 30},
      {[0, 0, -5, -6], -40}
    ],
    column_order = [1, 2, 3, 4]
  },
  {ok, V2} = gauss_solution(A2),

  [
    ?_assertEqual({error, no_simple_solution}, gauss_solution(A0)),
    ?_assertEqual(true, Check_solution(A1, V1)),
    ?_assertEqual(true, Check_solution(A2, V2))
  ].

get_sol_test_() ->
  X = [A * 0.0001 || A <- lists:seq(1, 50000)],

  V1 = [0, 0, 1, 0],
  V2 = [0.5, 2.5, 4.5, 0],
  V3 = [1, -0.5, -2.5, 0, 10, 100],

  Base = [V1, V2, V3],
  Map_fun =
    fun(V) ->
      P = [{X_val, tapol_epol:calc_val(V, X_val)} || X_val <- X],
      {ok, S} = get_least_squares_solution(P, length(V) - 1),
      ?_assertEqual(true, lists:all(fun({A, B}) -> abs(A - B) < 0.01 end, lists:zip(V, S)))
      end,
  lists:map(Map_fun, Base).

-endif.