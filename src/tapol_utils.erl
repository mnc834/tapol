%%%-------------------------------------------------------------------
%%% @author PKQ874
%%% @copyright (C) 2016, <COMPANY>
%%% @doc
%%%
%%% @end
%%% Created : 25. апр 2016 15:17
%%% module for keeping aux functions
%%%-------------------------------------------------------------------
-module(tapol_utils).
-author("PKQ874").

-include("tapol.hrl").

%% API
-export([are_equal/2]).

-spec are_equal(A :: float(), B :: float()) -> boolean().
%% checks if the absolute diferense between the wo floats is less than ?EPSOLINE
are_equal(A, B) ->
  abs(A - B) < ?EPSILON.
