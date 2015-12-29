(* -*-OCaml-*-
********************************************************************************
*
* File:         tree.ml
* RCS:          $Header: $
* Description:  
* Author:       Staal Vinterbo
* Created:      Sat Nov  5 17:27:44 2011
* Modified:     Sat Oct 10 18:30:51 2015 (Staal Vinterbo) staal@klump.gateway.pace.com
* Language:     caml
* Package:      N/A
* Status:       Experimental
*
* tree.ml is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* tree.ml is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with tree.ml; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*
* (c) Copyright 2011, Staal Vinterbo, all rights reserved.
*
********************************************************************************
*
* Revisions:
*
* Sat Oct 10 18:30:37 2015 (Staal Vinterbo) staal@klump.gateway.pace.com
*  Batteries syntax changed: -| to %
********************************************************************************
*)

open Batteries;;
open Genlex;;
open Utils;;


(* --------- Tree-related code ------------ *)

type tree = Node of (int * tree list);;

(* tree fold *)
let rec treefold f g z tree =
  let rec treefolds f g z subtrees =
    match subtrees with
      [] -> z
    | (x::xs) -> g (treefold f g z x) (treefolds f g z xs) in
  match tree with Node(label, subtrees) -> f label (treefolds f g z subtrees);;

let preorder = treefold List.cons (@) [];;

let edgef label ll =
  let f (i, pl) res = (label, i) :: (pl @ res) in
  (label, List.fold_right f ll []);;

let preedges tree = let (_, edges) = treefold edgef List.cons [] tree in edges;;

(* extract root of tree *)
let root tree = match tree with Node(r, _) -> r;;

(* # vertices in tree *)
let order = treefold (fun _ x -> x + 1) (+) 0;;

(* parse representation of tree (using Genlex and camplp4) *)
let lexer = make_lexer ["("; ")"; ","; "["; "]"];;
let rec parse_tree = parser
  | [< 'Kwd "("; 'Int label; 'Kwd ","; 'Kwd "["; nodes=parse_nodelist; 'Kwd "]" ;'Kwd ")" >] ->
      Node(label, nodes)
and parse_nodelist = parser
  | [< node = parse_tree ; rest = parse_tail >] -> node::rest
  | [< >] -> []
and parse_tail = parser
  | [< 'Kwd "," ; node = parse_tree; rest = parse_tail >] -> node::rest
  | [< >] -> [];;
let tree_from_string = parse_tree % lexer % Stream.of_string;;

(** color coding stack operations: lets us access next edge parent in constant time *)

type operations = PUSH | POP | POPPUSH | C | TOP;;

(* implementations of operations *)
let opfun x =
  match x with
    PUSH -> (fun l current -> (current::l, current))
  | POPPUSH -> (fun (x::xs) current -> (current::xs, current))
  | POP -> (fun (x::y::ys) _ -> (y::ys, y))
  | C -> (fun l current -> (l, current))
  | TOP -> (fun (x::xs) current -> (x::xs, x));;

(* fix rightmost subtree: repair (using splitwhile). *) 
(*   Translates first PUSH and TOP into POPPUSH and POP *)
let rec splitwhile p l =
  match l with
    | [] -> ([], [])
    | (x::xs) ->
      if p x then let (left, right) = splitwhile p xs in (x::left, right)
      else ([], l);;

let repair oplist =
  let (left, right) = splitwhile (fun x -> x == C) oplist in
  match right with
    | [] -> oplist
    | (x::xs) -> match x with
	| PUSH -> left @ (POPPUSH::xs)
        | TOP -> left @ (POP::xs)
        | _ -> oplist;;


(* treefold f function for generating stack operations *)
let opsf label ilist =
  match ilist with
    |  [] -> [TOP]
    | (l::[]) -> C::l
    | _ -> let (x::xs) = List.rev ilist in
           PUSH::(List.concat (List.rev ((repair x) :: xs)));;

let treeopsl = treefold opsf List.cons [];;
let treeops tree = List.map opfun (treeopsl tree);;

(* initial stack *)
let initstack v = [v;v];;

(* preorder (operation,edge) list from tree *)
let ielist tree = let (x::xs) = treeops tree in zip xs (preedges tree);;
