(* -*-OCaml-*-
********************************************************************************
*
* File:         yagmproto.ml
* RCS:          $Header: $
* Description:  An Ocaml implementation of matching a list of trees in a graph
*               using color coding. By graphs we mean simple attributed graphs.
*               Input is provided in three files, the first containing the trees, 
*               each tree is given on a single line in a recursive
*               "(vertex index, [subtrees])"  format.
*               The next two files contain real valued matrices 
*               A and W in comma separated row major layout. The A matrix
*               is a matrix where entry i,j contains the similarity of 
*               vertex i in a graph G1 with vertex j in Graph G2. Graph G1
*               must be a superset of all the trees given. Graph G2 is the
*               target graph in which we seek to match the trees.  
*               The matrix W is the edge similarity matrix  
*               of G1 and G2. Edges for each graph are enumerated
*               as in the row major linear storage of the adjacency 
*               matrix of the graph. In any case, self loops, i.e., 
*               the diagonals are not counted. In the undirected case,
*               when the adjacancy matrix is symmetric, 
*               the upper triangular entries are also not counted.
*               Alternatively the matrices A and W contain the graphs
*               G1 and G2 respectively. The first row contains the weights
*               associated with vertices, the subsequent rows are the weighted
*               adjacency matrix of the graph.
*
* Author:       Staal Vinterbo
* Created:      Thu Aug  4 10:12:08 2011
* Modified:     Mon Oct 31 20:23:54 2011 (Staal Vinterbo) staal@dink
* Language:     caml
* Package:      N/A
* Status:       Experimental
*
* yagmproto.ml is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* yagmproto.ml is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with yagmproto.ml; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*
* (c) Copyright 2011, Staal Vinterbo, all rights reserved.
*
********************************************************************************
*)

open Batteries;;
open Utils;;
open Tree;;
open State;;
open Ga;;

(********** Color coding as a fold over preorder tree edges *)

(* a candidate tree consists of ((endpoint, colorset), (score, vertices in reverse, endpoint stack) *)

(* extend a candidate tree ptree using origin exte edge *)
let extend nextlen state bound edge ptree =
  let (op, (a, b)), ((colorset,u), (score, vlist, stack)) = edge, ptree in
  let extrip v =
    if (colorset land state.color.(v)) == 0 then
      let q = score +. (state.v b v) +. (state.w (a,b) (u,v)) in
      if state.keep nextlen bound q then
	let (ns, nu) = op stack v in
	Some ((colorset lor state.color.(v), nu),(q, v::vlist, ns))
      else None
    else None in
  Enum.filter_map extrip (Array.enum state.tadj.(u));;

(* fold operator for single trial which is a fold over a list of preorder edges with stack instructions *)
let oplus bound (nextlen, state, table) edge =
  let (n,m,k) = (state.n, state.m, nextlen - 1) in
  let op table ptree = Enum.fold state.store.update table (extend nextlen state bound edge ptree) in
  (nextlen + 1, state, Enum.fold op (state.store.refresh (n,m,nextlen) table) (state.store.items (n,m,k) table));;

(* extract a list of tree vertex lists that have opt score *) 
let emax ?compare:(compare=Pervasives.compare) e =
  Enum.fold (fun (qq,ll) (q, l, _) ->
    match compare qq q with
      -1 -> (q,[l])
    | 0 -> (qq, l::ll)
    | _ -> (qq, ll)) (0., []) e;;

(* single trial: fold oplus over the (instruction,edge) list, return best matches *)
let singletrial state bound ilist = 
  let init = [? ((state.color.(u),u), (state.v state.root u, [u], initstack u)) |
                u <- List: state.gvertices  ?] in
  let initf = Enum.filter (fun (key, (q, l, s)) -> state.keep 1 bound q) init in
  let table = Enum.fold state.store.update (state.store.create (state.n, state.m, 1) state.store) initf in
  let _,_,ftable = List.fold_left (oplus bound) (2, state, table) ilist in
  emax ~compare:state.compare (state.store.values (state.n, state.m, state.m) ftable);;

(* treematch: perform predetermined number of trials, returning the highest score match *)
let treematch ?bound:(bound=0.0) state tree =
  let trial ilist (boundt,ll) _ =
    let () = state.color <- state.new_colors () in
    let (q, l) = singletrial state boundt ilist in
    match state.compare boundt q with
      -1 -> (q,l)
    | 0 -> (boundt, l @ ll)
    | _ -> (boundt, ll)
  in
  let ilist = ielist tree in
  let (score, l) = Enum.fold (trial ilist) (bound, []) (0 --^ state.times) in
  (score, List.map (state.make_assignment % List.rev) l);;


(* Yagma algorithm ********************* *)

(* -------- GA Fitness Function ------- *)
module IntMap = Map.Make(struct type t = int let compare = compare end);;
  
let ascore v w tedges tnodes ind =
  try 
    let l = [? List : (a,b) | (a,b) <- List : (List.map2 (fun x y -> (x,y)) tnodes (Array.to_list ind)) ?] in
    let lmap = IntMap.of_enum [? (a,b) | (a,b) <- List: l ?] in
    let fl i = IntMap.find i lmap in
    let ascore = List.fold_left (+.) 0.0 (List.map (uncurry v) l) and
	wscore = List.fold_left (+.) 0.0 (List.map (fun (a,b) -> w (a, fl a) (b, fl b)) tedges) in
    ascore +. wscore
  with Not_found -> 0.0;;
let flipw w (a1,b1) (a2,b2) = w (a1, a2) (b1, b2);;
let ga_bound verbose n state tree =
      let tedges = preedges tree and
	  tnodes = preorder tree in
      let gafit = ascore state.v (flipw state.w) tedges tnodes in
      let () = if verbose then debug " (running GA..." and
	  gabound, _ = runaga 1000 128 (makeri gafit (order tree) n) gafit in
      let  () = if verbose then debug ("bound: " ^ (Printf.sprintf "%f" gabound) ^ ") ") in
      gabound;;


(* ----- Create weighted bi-partite graph of tree match assignments ------ *)
type inputmode = Similarity | Adjacency;;
  
let graphmatch ?verbose:(verbose=false) ?gbound:(gbound=false) ?mode:(mode=Similarity) ?strict:(strict=true) ?bound:(bound=0.) eps amatrix wmatrix trees =
  let store = hash_store in
  let m, n =
    (if mode == Similarity then Array.length amatrix, Array.length amatrix.(0) else
    Array.length amatrix.(0), Array.length wmatrix.(0)) in
  let bigraph = Array.make_matrix m n 0.0 in
  let dotree i tree =
    let state =
      (if mode == Similarity then
        create_state_sim ~strict:strict store
      else create_state_adj ~strict:strict store) eps wmatrix amatrix tree in
    let start_bound = if gbound then
          (ga_bound verbose n state tree) -. 1.0e-10
    else
      bound in
    let (score, matches) = treematch ~bound:start_bound state tree in
    begin
      List.iter (fun (a,b) -> bigraph.(a).(b) <- bigraph.(a).(b) +. score) (List.concat matches);
      if verbose then begin
	let mlen = List.length matches in
	debug ("tree # " ^ (Int.to_string i) ^ "  score: " ^ (Float.to_string score) ^
	       ", # matches: " ^ (Int.to_string mlen) ^ "\n");
      end;
    end;
  in
  Enum.iteri dotree trees;
  bigraph;;


(* -------------------------- Program entry point ------------------------ *)
					      
let version = "1.36";;
let main() = 
  flush stderr;
  let eps = ref 0.2 and
      verbose = ref false and
      seed = ref (int_of_float (Unix.time ())) and
      gbound = ref false and
      bound = ref 0. and
      strict = ref false and
      mode = ref Similarity and
      file1 = ref "" and
      file2 = ref "" and
      file3 = ref "" and
      files = ref [] and
      usage = (Sys.argv.(0) ^ " version " ^ version ^
                 ". (c) 2011-2015 Staal A. Vinterbo.\n" ^
	       "usage: " ^ Sys.argv.(0) ^
               " [-e FLOAT] [-s INT] [-b FLOAT] [-v] [-g] TFILE AFILE WFILE") in
  let speclist = [ ("-e", Arg.Float  (fun ein -> eps := ein),
                    ": set eps to FLOAT.");
		   ("-b", Arg.Float  (fun ein -> bound := ein),
                    ": set the starting lower bound to FLOAT.");
                   ("-g", Arg.Unit    (fun () -> gbound := true;()),
                    ": compute initial bounds using GA.");
                   ("-strict", Arg.Unit    (fun () -> strict := true;()),
                    ": stricter pruning of search space.");
                   ("-m", Arg.Unit    (fun () -> mode := Adjacency;()),
                    ": files contain graph matrices (first row is vertex weights,\n"^
                      "       while following rows contain the weighted adjacency matrix).");
		   ("-seed", Arg.Int    (fun sin -> seed := sin),
                    ": seed random number generator with INT.");
		   ("-v", Arg.Unit    (fun () -> verbose := true;()),
                    ": print progress notes.");
		   ("-stdin", Arg.Unit (fun () -> files := ("-"::!files);()),
		    ": specify stdin as input for this positional file argument.") 
		 ] in
  let () = begin
    Arg.parse speclist (fun x -> files := (x::!files)) usage;
    if List.length !files != 3 then begin
      Printf.eprintf "Need three filenames. Try -help for usage information.\n";
      raise (Arg.Bad "No filenames given.");
    end;
    Random.init !seed;
    file1 := List.last !files;
    file2 := List.nth !files 1;
    file3 := List.first !files;
    if !verbose then
      debug ("Reading data from " ^ !file1 ^ " " ^ !file2 ^ " " ^ !file3 ^
	     ".\nParameters: eps = " ^ (Float.to_string !eps) ^
             ", seed = " ^ (Int.to_string !seed) ^
             ", bound = " ^ (Float.to_string !bound) ^
             ", strict = " ^ (if !strict then "yes" else "no") ^
             ", input mode = " ^ (if !mode == Similarity
               then "Similarity" else "Adjacency") ^ ".\n");
  end in
  let trees = Enum.map tree_from_string (readlines !file1) and
      amatrix = readmatrix !file2 and
      wmatrix = readmatrix !file3 in
  begin
    if !verbose then debug("Starting computations, please wait..\n");
    let stime = Unix.time () in 
    let result = graphmatch ~verbose:!verbose ~gbound:!gbound ~mode:!mode ~strict:!strict ~bound:!bound !eps amatrix wmatrix trees in begin
      if !verbose then
        debug("Computations took: " ^
                 (Float.to_string ((Unix.time()) -. stime)) ^ " seconds.\n");
      printmatrix "%f" result;
    end
  end;;


(* --- call the entry point --- *)
main();;


