(* -*-OCaml-*-
********************************************************************************
*
* File:         state.ml
* RCS:          $Header: $
* Description:  
* Author:       Staal Vinterbo
* Created:      Sat Nov  5 17:37:37 2011
* Modified:     Sat Oct 10 19:22:43 2015 (Staal Vinterbo) staal@klump.gateway.pace.com
* Language:     caml
* Package:      N/A
* Status:       Experimental
*
* state.ml is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* state.ml is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with state.ml; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*
* (c) Copyright 2011, Staal Vinterbo, all rights reserved.
*
********************************************************************************
*
* Revisions:
*
* Sat Oct 10 19:21:36 2015 (Staal Vinterbo) staal@klump.gateway.pace.com
*  updated to batteries 2, List.sort's ~cmp is no longer optional, and 
*  -| is now %.
********************************************************************************
*)

open Batteries;;
open Utils;;
open Tree;;



(* abstracted store access ********** *)
(* needs to define:
   update table (key, value) -> table
   items (n,m,k) table -> items in table with colorsets of size k (out of m)
   values (n,m,k) table -> (key, values) in table with colorsets of size k (out of m)
   refresh (n,m,k) table -> "refreshed table"
   create (n,m,k) -> new table
   size (n,m,k) table -> the number of elements in the table with colorsets of size k (out of m)
   estimate_size (n,m,k) -> estimate the size of needed storage
*)

type 'a store = {
  update : 'a -> ((int * int) * (float * int list * int list)) -> 'a ;
  items : (int * int * int) -> 'a -> ((int * int) * (float * int list * int list)) Enum.t ;
  values: (int * int * int) -> 'a -> (float * int list * int list) Enum.t ;
  estimate_size: (int * int * int) -> int ;
  refresh: (int * int * int) -> 'a -> 'a ;
  create: (int * int * int) -> 'a store -> 'a ;
  size: (int * int * int) -> 'a -> int;
  mutable table: 'a
};;


(* estimate size of storage needed for colorsets of length len and G order n *)

(* binomial coefficient *)
let binom n k =
  let k' = min (n - k) k and
      r = ref 1 and
      n' = ref n in
  for d = 1 to k' do
    r := !r * !n'; n' := !n' - 1; r := !r / d;
  done;
  !r;;

let max_size n m k = n * (binom m k);;

(* provide update, items, values and refresh operations on some map data store.
   here Hashtbl *)

let h_update h ((s, u),(q, l, ss)) =
  match Hashtbl.find_option h (s, u) with
    None -> begin Hashtbl.add h (s, u) (q, l, ss); h end
  | Some (q',l', s') -> begin if (q > q') then Hashtbl.replace h (s,u) (q,l,ss); h end;;
let h_items (n,m,k) = Hashtbl.enum;; (* all key,value pairs *)
let h_values x = Enum.map snd % (h_items x) ;; (* all the values, without keys *)
let h_estimate_size (n,m,k) = (max_size 1 m k);;
let h_create (n,m,k) store =  Hashtbl.create (h_estimate_size (n,m,k));;
let h_refresh (n,m,k) table = Hashtbl.create (h_estimate_size (n,m,k));; (* create a new store *)
let h_size x = Hashtbl.length;;



(* Large static table, allocated once *)
(* given a colorset of size k, compute the next size k colorset in
   the lexicographic order of all size k out of n colorsets *)
let next_c x =
  let u = x land (-x) in
  let v = u + x in
  if (v == 0) then 0 (* overflow *) else v + (((v lxor x)/u) lsr 2);;

(* generate enum of all colorsets of size k out of n colors *)
let colorsets n k =
  Enum.fold (fun (latest::l) _ ->
    (next_c latest)::latest::l) [(1 lsl k) - 1] (1 --^ (binom n k));;

let t_update h ((s, u),(q, l, ss)) =
  (*debug("update " ^ (Int.to_string s) ^ ", " ^ (Int.to_string u) ^ "\n"); *)
  let (q', _, _) = h.(s-1).(u) in
  if q > q' then
    h.(s-1).(u) <- (q,l,ss);
  h;;

let t_items (n,m,k) h =
  Enum.filter_map
    (fun (key, (q,l,ss)) -> match q with
        0. -> None
      | _ -> Some (key, (q,l,ss)))
    [? ((s,u), h.(s-1).(u)) | s <- List : colorsets m k; u <- 0 --^ n ?];;
let t_values x = Enum.map snd % (t_items x) ;;
let t_estimate_size (n,m,k) = (max_size n m k);;
let t_refresh (n,m,k) table = table;;
let t_create (n, m, k) store =
  let rows = (1 lsl m) - 1 and
      cols = n in Array.make_matrix rows cols (0., [], []);; (*
  if (Array.length store.table) == rows && (Array.length store.table.(0)) == cols then
    begin
      Array.iter (fun a -> Array.fill a 0 cols (0., [], [])) store.table;
      store.table
    end
  else let tab = Array.make_matrix rows cols (0., [], []) in
       begin
         store.table <- tab;
         tab
       end;;
                                                             *)
let t_size x m = match m with [||] -> 0 | _ -> (Array.length m) * (Array.length m.(0));;


let hash_store = {
  update = h_update;
  items = h_items;
  values = h_values;
  estimate_size = h_estimate_size;
  refresh = h_refresh;
  create = h_create;
  size = h_size;
  table = Hashtbl.create 0;
};;

let array_store = {
  update = t_update;
  items = t_items;
  values = t_values;
  estimate_size = t_estimate_size;
  refresh = t_refresh;
  create = t_create;
  size = t_size;
  table = Array.make_matrix 1 1 (0., [], [])
};;


(*  Color coding as a fold over tree edges in preorder ************* *)
(*  :  foldl oplus table edges *)

(* state threaded through calls (monads anyone?)  *)
type 'a state = {
    times : int;
    gvertices : int list;
    root : int;
    m : int;
    n : int;
    v : (int -> int -> float);
    w : ((int*int) -> (int*int) -> float);
    new_colors : (unit -> (int array));
    mutable color : int array;
    keep : (int -> float -> float -> bool);
    make_assignment : ((int list) -> (int * int) list);
    tadj : int array array;
    compare: float -> float -> int;
    store: 'a store
  } ;;

(* ------------ Messy computing of 'state' -------- *)

(*  The inputs are two matrices amatrix and wmatrix.
    wmatrix is an edge similarity matrix, i.e, 
    wmatrix.(i).(j) contains the similarity between edge i
    in graph 1 and edge j in graph 2. Edges are enumerated
    as in the row major linear storage of the adjacency 
    matrix of the graph. In any case, self loops, i.e., 
    the diagonals are not counted. In the undirected case,
    when the adjacancy matrix is symmetric, 
    the upper triangular entries are also not counted.
    Whether graphs are directed or not is determined by 
    the dimensions of wmatrix which correspond to the number of 
    counted edges (in the complete simple graphs).
    The parameters m and n correspond to the number of vertices
    in graph 1 and graph 2, respectively.
*)

(* figure out whether the graphs are directed or not *)
let directedp m n wmatrix =
  let nrow,ncol = Array.length wmatrix, Array.length wmatrix.(0) in
  (nrow == (intexp m 2) - m,  ncol == (intexp n 2) - n);;

(* map edge (i,j) to index in wmatrix *)
let eidx m n wmatrix =
  let (d1, d2) = directedp m n wmatrix and
      f1d = (fun i j -> let c = i*m + j in c - (c/(m + 1) + 1)) and
      f2d = (fun i j -> let c = i*n + j in c - (c/(n + 1) + 1)) and
      fu = (fun i j ->
	let (x,y) = (if i > j then (i,j) else (j,i)) in (x*(x - 1))/2 + y) in
  ((if d1 then f1d else fu), (if d2 then f2d else fu));;

let make_v amatrix = (fun a b -> amatrix.(a).(b));;

let make_w m n wmatrix =
  let edgei1, edgei2 = eidx m n wmatrix in
  (fun (a1, a2) (b1, b2) ->
    wmatrix.(edgei1 a1 a2).(edgei2 b1 b2));;

  
(* -- Create initial colorsets --*)
let make_colors n k =
  [? Array : 1 lsl x | x <- Enum.map (fun _ -> Random.choice (0 --^ k)) (0 --^ n) ?];;
    

(* --- Pruning of solutions --- *) 
(* -- compute array that estimates the max additional value an assignment can
      get when extended by i tuples --*)
let boundsarray amatrix wmatrix tree =
  let m, n = Array.length amatrix, Array.length amatrix.(0) and
      tnodes = preorder tree and
      tpairs = preedges tree in
  let ei1, ei2 = eidx m n wmatrix and
      k = List.length tnodes in
  let teidxs = List.map (uncurry ei1) tpairs in
  let widxs = product teidxs [? List : x | x <- 0 --^ (Array.length wmatrix.(0)) ?] and
      aidxs = product tnodes [? List : x | x <- 0 --^ n ?] in
  let wvals = List.take k (List.sort (fun x y -> compare y x) (List.map (fun (i,j) -> wmatrix.(i).(j)) widxs)) and
      avals = List.take k (List.sort (fun x y -> compare y x) (List.map (fun (i,j) -> amatrix.(i).(j)) aidxs)) in
  let tmp = addflp wvals avals in
  Array.rev (cumsum (+.) 0. (0.::tmp));;

(* -- function to create function for filtering sub-matches  -- *)
let keep_p ?strict:(strict=true) barray nextlen bound q =
  if strict then bound < q +. barray.(nextlen) else bound <= q +. barray.(nextlen);;

(* - list of edges in an undirected complete graph of order n *)
let kedgesudir n =
  let l = List.of_enum (0 --^ n) in
  List.rev [? List : (x,y) | (x,y) <- List : product l l; x > y ?];;

(* - list of edges in a directed complete graph of order n *)
let kedgesdir n =
  let l = List.of_enum (0 --^ n) in
  List.rev [? List : (x,y) | (x,y) <- List : product l l; x != y ?];;

(* - indices of elements in array a for which p holds *)
let which p a = filter (fun i -> p a.(i)) (0 --^ (Array.length a));;

(* extract an adjacency matrix (array of arrays representation) from watrix *)
let adjlist sel m n wmatrix =
  let dt, dg = directedp m n wmatrix in
  let eit, eig = eidx m n wmatrix in
  let row m i = m.(i) and
      col m i = Array.map (fun r -> r.(i)) m in
  let k,ei,d,acc = sel ((m,eit,dt,row), (n,eig,dg,col)) in
  let edges = (if d then kedgesdir else kedgesdir) k and
      good a = (Array.fold_left (+.) 0. a) > 0.0 and
      am = Array.make_matrix k k 0 in
  begin
    List.iter (fun (a,b) ->
      let i = ei a b in
	am.(a).(b) <- (if (good (acc wmatrix i)) then 1 else 0)) edges;
    (*printmatrix "%d" am;*)
    Array.map (fun row -> Array.of_enum (which (fun x -> x > 0) row)) am;
  end;;

let oadjlist = adjlist fst;;  (* adjlist for origin *)
let tadjlist = adjlist snd;;  (* adjlist for target *)


let sim a b = 1. -. (Float.abs (a -. b)) (*
  if a == b then 1.
  else 1. -. (Float.abs (a -. b)) /. (max (Float.abs(a)) (Float.abs(b))) *);;

let make_v2 amatrix wmatrix =
  (fun t u -> sim amatrix.(0).(t) wmatrix.(0).(u));;

let make_w2 amatrix wmatrix =
  (fun (t,t') (u, u') -> sim amatrix.(t + 1).(t') wmatrix.(u + 1).(u'));;

let flip f x y = f y x;;
let boundsarray2 gedges tedges gnodes tnodes v w =
  let k = List.length tnodes in
  let taketop l = List.take k (List.sort (flip compare) l) in 
  let avals =
    taketop [? List: v t u | t <- List: tnodes; u <- List: gnodes ?] and
      wvals =
    taketop [? List: w a b | a <- List: tedges; b <- List: gedges ?] in
  let tmp = addflp wvals avals in
  let ba = Array.rev (cumsum (+.) 0. (0.::tmp)) in
  begin
    (*debug ("barray " ^
             (Enum.fold (fun s x -> s ^ " " ^ (Float.to_string x)) ""
                (Array.enum ba)) ^ "\n"); *)
    ba
  end;;

let adjlist2 gmat =
  let n = Array.length gmat.(0) in
  let nbors arr =
      Array.of_enum (Enum.filter_map (fun i ->
        match arr.(i) with 0. -> None | _ -> Some i) (0--^n)) in
  Array.of_enum (Enum.map (fun i -> nbors gmat.(i)) (1--n));;

let graphedges admat =
  let edgs i a = [? List : (i, x) | x <- Array: a ?] in
  Enum.foldi (fun i a l -> (edgs i a) @ l) [] (Array.enum admat);;


(* compute the 'state' for a given program input and tree -- *)
let create_state_sim ?strict:(strict=true) store eps wmatrix amatrix tree =
  let m = Array.length amatrix and
      n = Array.length amatrix.(0) and
      barray = boundsarray amatrix wmatrix tree and
      k = order tree in
  let tvertices = preorder tree in {
    times = int_of_float ((log (1.0 /. eps)) *. (exp (float_of_int k)));
    gvertices = [? List : x | x <- 0 --^ n ?];
    new_colors = (fun () -> make_colors n k);
    n = n;
    m = k;
    root = root tree;
    color = make_colors n k;
    v = make_v amatrix;
    w =  make_w m n wmatrix;
    keep  = keep_p ~strict:strict barray; 
    make_assignment = (fun l -> zip tvertices l);
    tadj = tadjlist m n wmatrix;
    compare = Pervasives.compare;
    store = store;
  };;

let create_state_adj ?strict:(strict=true) store eps wmatrix amatrix tree =
  let n = Array.length wmatrix.(0) and
      v = make_v2 amatrix wmatrix and
      w = make_w2 amatrix wmatrix and
      tadj = adjlist2 wmatrix and
      tedges = preedges tree and
      tnodes = preorder tree in
  let gnodes = [? List : x | x <- 0 --^ n ?] and
      gedges = List.filter (fun (a,b) -> a != b) (graphedges tadj) in
  let barray = boundsarray2 gedges tedges gnodes tnodes v w and
      k = order tree in {
        times = int_of_float ((log (1.0 /. eps)) *. (exp (float_of_int k)));
        gvertices = [? List : x | x <- 0 --^ n ?];
        new_colors = (fun () -> make_colors n k);
        n = n;
        m = k;
        root = root tree;
        color = make_colors n k;
        v = v;
        w = w;
        keep  = keep_p ~strict:strict barray; 
        make_assignment = (fun l -> zip tnodes l);
        tadj = tadj;
        compare = Pervasives.compare;
        store = store;
      };;
