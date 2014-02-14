(* -*-OCaml-*-
********************************************************************************
*
* File:         utils.ml
* RCS:          $Header: $
* Description:  yagm utilities definitions
* Author:       Staal Vinterbo
* Created:      Sat Nov  5 16:54:56 2011
* Modified:     Sun Nov  6 16:24:14 2011 (Staal Vinterbo) staal@mats
* Language:     caml
* Package:      N/A
* Status:       Experimental
*
* utils.ml is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* utils.ml is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with utils.ml; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*
* (c) Copyright 2011, Staal Vinterbo, all rights reserved.
*
********************************************************************************
*)


open Batteries;;
(* ------------------ Utility definitions --------------------- *)

(* List.cartesian_product *)
let product xs ys =
  List.fold_left (fun acc x ->
    List.fold_left (fun acc y -> (x,y) :: acc) acc ys) [] xs;;

(* cumulative "sum" of elements in a list, returned as an array *)
let cumsum op zero list =
  Array.of_list (List.rev (snd (List.fold_left
				  (fun (x,l) y ->
				    let res = op x y in
				    (res, res::l)) (zero,[]) list)));;

let uncurry f (a,b) = f a b;;

(* debug print function *)
let debug s = prerr_string s;flush stderr;;

(* exponentiations in ocaml is only for floats *)
let intexp a b = int_of_float ((float_of_int a) ** (float_of_int b));;

let zip a b =
  let rec ziph acc a b =
    match a,b with
      [], _ -> acc
    | _, [] -> acc
    | (x::xs),(y::ys) -> ziph ((x,y)::acc) xs ys in
  List.rev (ziph [] a b);;

let addtwo op z a b =
  let rec a2 acc a b =
    match a,b with
      [],[] -> acc
    | (x::xs),[] -> a2 ((op x z)::acc) xs b
    | [], (y::ys) -> a2 ((op z y)::acc) a ys
    | (x::xs), (y::ys) -> a2 ((op x y)::acc) xs ys in
  List.rev(a2 [] a b);;

let addflp = addtwo (+.) 0.;;


(* IO ----------- *)

(* read an enumeration of lines from a file, '-' as filename means stdin *)
let readlines ?p:(p=(fun s -> (String.length s) > 0 && s.[0] != '#')) file =
  let input = match file with
    "-" -> stdin
  | _ -> open_in file in
  Enum.filter p (IO.lines_of input);;

(* read a comma separated matrix of floats from file with name 'file' *)
let readmatrix file =
  let s2a s = [? Array : (String.to_float (String.trim ss)) |
          ss <- List : (String.nsplit s ",") ?] in
  try
    [? Array : (s2a s) | s <- (readlines file) ?]
  with _ -> raise (Failure file);;


(* print a comma separated matrix *)
let fprint fmt stream x = IO.printf stream fmt x;;
let printmatrix ?stream:(stream=IO.stdout) fmt matrix =
  Array.iter (fun a ->
    Array.print ~first:"" ~sep:"," ~last:"\n" (fprint fmt) stream a) matrix;;

