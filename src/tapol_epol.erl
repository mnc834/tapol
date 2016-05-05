%%%-------------------------------------------------------------------
%%% @author PKQ874
%%% @copyright (C) 2016, <COMPANY>
%%% @doc
%%%
%%% @end
%%% Created : 25. апр 2016 14:58
%%% a library for working with polynomials
%%%-------------------------------------------------------------------
-module(tapol_epol).
-author("PKQ874").

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

-type e_polynomial() :: [float()].
%% @doc Polynomial coefficients, the higherst's power coefficient comes first:
%% P(X) = An * X^n + An-1 * X^(n-1) + .. + A0 yields coefficients:
%% [An, An-1, .. A0]

-export_type([e_polynomial/0]).

%% API
-export([calc_val/2,
  derivative/1,
  stretch/2,
  add/2
]).

-spec calc_val(P :: e_polynomial(), X :: float()) -> float().
%% @doc calculates value of polynomial P in the point X,
%%      for an empty polynomials throughs exception empty_polynomial
calc_val([], _) ->
  erlang:error(empty_polynomial);
calc_val([H | T], X) ->
  Foldl_fun =
    fun(A, S) ->
      S * X + A
    end,
  lists:foldl(Foldl_fun, H, T).

-spec derivative(P :: e_polynomial()) -> e_polynomial().
%% @doc calculates a derivative for the given polynomial P
%%      for an empty polynomials throughs exception empty_polynomial
derivative([]) ->
  erlang:error(empty_polynomial);
derivative(P) ->
  Foldr_fun =
    fun(C, N) ->
      {C * N, N + 1}
    end,
  {L, N} = lists:mapfoldr(Foldr_fun, 0, P),
  {D_p, _} = lists:split(N - 1, L),
  D_p.

-spec stretch(P :: e_polynomial(), S :: float()) -> e_polynomial().
%% @doc implements streching by substituting the arguments X = S * X
%%      returns the new coefficients
stretch(P, S) ->
  Foldr_fun =
    fun(C, V) ->
      {C * V, V * S}
    end,
  {L, _} = lists:mapfoldr(Foldr_fun, 1, P),
  L.

-spec add(P1 :: e_polynomial(), P2 :: e_polynomial()) -> e_polynomial().
%% @doc adds polynomials
add(P1, P2) ->
  N1 = length(P1),
  N2 = length(P2),

  {L1, L2} =
    case N1 =:= N2 of
      true ->
        {P1, P2};
      false ->
        N = abs(N1 - N2),
        L = [0 || _ <- lists:seq(1, N)],
        case N1 > N2 of
          true ->
            {P1, lists:append(L, P2)};
          false ->
            {lists:append(L, P1), P2}
        end
    end,
  lists:map(fun({A, B}) -> A + B end, lists:zip(L1, L2)).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Unit tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-ifdef(TEST).

are_polynomials_equal(P1, P2) when is_list(P1), is_list(P2) ->
  lists:all(fun({A, B}) -> tapol_utils:are_equal(A, B) end, lists:zip(P1, P2));
are_polynomials_equal(A1, A2) ->
  A1 =:= A2.

calc_val_test_() ->
  X = 100,

  Res_1 =
    try calc_val([], X) of
      Val_1 ->
        Val_1
    catch
      error : Error ->
        Error
    end,
  Res_2 = calc_val([0], X),
  Res_3 = calc_val([1, 1, 1], X),
  [
    ?_assertEqual(empty_polynomial, Res_1),
    ?_assertEqual(true, tapol_utils:are_equal(0, Res_2)),
    ?_assertEqual(true, tapol_utils:are_equal(X * X + X + 1, Res_3))
  ].

derivative_test_() ->
  P = [10.4, 1, 2, 3, 4],

  Derivatives =
    [
      [41.6, 3, 4, 3],
      [124.8, 6, 4],
      [249.6, 6],
      [249.6],
      [],
      empty_polynomial
    ],


  Res_fun =
    fun Derivative(P) ->
      try derivative(P) of
        D_p ->
          [D_p | Derivative(D_p)]
      catch
        error : Error ->
          [Error]
      end
    end,
  Result = Res_fun(P),

  [?_assertEqual(true, are_polynomials_equal(P1, P2)) || {P1, P2} <- lists:zip(Derivatives, Result)].

strech_test_() ->
  P = [1, 2, 3],
  S = 10.0,
  X = 11.1,
  P_s = [100, 20, 3],

  Res_1 = stretch(P, S),
  Res_2 = calc_val(P, S * X),
  Res_3 = calc_val(Res_1, X),
  [
    ?_assertEqual(true, are_polynomials_equal(Res_1, P_s)),
    ?_assertEqual(true, tapol_utils:are_equal(Res_2, Res_3))
  ].

add_test_() ->
  P1 = [1, 2, 3, 4, 5],
  P2 = [4.1, 5.2, -6.3],
  Res = [1, 2, 7.1, 9.2, -1.3],

  Res_1 = add(P1, P2),
  Res_2 = add(P2, P1),
  Res_3 = add([], P1),
  Res_4 = add(P2, []),
  Res_5 = add([], []),

  [
    ?_assertEqual(true, are_polynomials_equal(Res, Res_1)),
    ?_assertEqual(true, are_polynomials_equal(Res, Res_2)),
    ?_assertEqual(true, are_polynomials_equal(P1, Res_3)),
    ?_assertEqual(true, are_polynomials_equal(P2, Res_4)),
    ?_assertEqual(true, are_polynomials_equal([], Res_5))
  ].

-endif.