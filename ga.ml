(* -*-OCaml-*-
********************************************************************************
*
* File:         ga.ml
* RCS:          $Header: $
* Description:  Genetic algorithm for the assignment problem
* Author:       Staal Vinterbo
* Created:      Sun Oct 30 19:58:02 2011
* Modified:     Sun Oct 30 19:59:10 2011 (Staal Vinterbo) staal@mats
* Language:     caml
* Package:      yagma
* Status:       Experimental
*
* ga.ml is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* ga.ml is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with ga.ml; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*
* (c) Copyright 2011, Staal Vinterbo, all rights reserved.
*
********************************************************************************
*)

open Batteries;;
open Std;;
(* ----------------- Ugly Genetic algorithm code --------------- *)

type indt = float * int array;;
type popt = indt list;;
type fitt = int array -> float;;

(* generate random permutation on n symbols *)
let randperm n =
  Random.shuffle (0--^n)

(* Array.take is missing from the standard library/Batteries *)
let atake k a = Array.sub a 0 k;;


(* generates a function that generates a random individual.
   Note that n must be smaller than 2^30. *)
let makeri (fitness:fitt) (m:int) (n:int) =
  (fun () -> let v = atake m (randperm n) in (fitness v, v));;


(* inverse of permutation *)
let invperm p = 
  let a = Array.make (Array.length p) (-1) in
  Array.iteri (fun i x -> a.(x) <- i) p;
  a;;

(* apply permutation *)
let permute what by = Array.map (fun where -> what.(where)) by;;

(* compute the preimage of an injection a as a Hashtbl.
   a represents an injection a(x) = a.(x). *)
let backedges a =
  let h = Hashtbl.create (Array.length a)
  in Array.iteri (fun i x -> Hashtbl.add h x i) a; h;;

(*  modified and extended cycle crossover for injections.

    we want to randomly mix two injections, creating an injection c
    such that
       c(x) in {a(x), b(x)}.

    Inputs are two lists of length m representing two injections
    a,b:{0,1,...,m} -> {0,1,...,n} for some n >= m (a sends i to a[i]).

    Algorithm sketch:
    The algorithm is an extension of cycle crossover with random selection
    of cycles (not altenating). This is extended to the injection case
    as explained below.

    If n = m we can view each of a and b as a perfect bipartite matching
    between Domain and Codomain.
    Let B be the bipartite graph constructed by the union of a and b.
    Each vertex in B has degree at most two, and B can be decomposed into
    disconnected components in which either all vertices have degree 1 or 2.
    All components where vertices have degree 1 are single edges, and
    represent the set of shared assigments. These are kept as they are.
    The remaining components are of vertice degree two, and are cycles.
    Of these cycles, we are only allowed to include every other edge for the
    result to be a perfect matching. Starting at a Domain vertex we can choose
    between the two directions we want to go, one direction includes
    edges from a, the other includes edges from b. We do this until we have
    traversed all cycles.

    If n > m, we can "ignore" cycle parts that have vertices outside
    of Domain and the union of Ranges. We simply stop processing an alternating
    path when we venture into the "ignored" part. The idea is to sample among
    all possible extensions of a and b into bijections and then processing
    as in the case where n = m.
*)
let mix a' b' =
  let perm = randperm (Array.length a') in
  let a = permute a' perm and
      b = permute b' perm in 
  let alen = Array.length a and
      used = Array.make (Array.length a) false and
      mixed = Array.make (Array.length a) (-1) and
      aback = backedges a and
      bback = backedges b and 
      start = ref 0 in
  while !start <> -1 & not used.(!start) do
    let (there,back) = if (Random.float 1.0) > 0.5 
      then (a, bback) else (b, aback) and
        x = ref !start
    in
    while !x <> -1 & not used.(!x) do
      mixed.(!x) <- there.(!x);
      used.(!x) <- true;
      x := if Hashtbl.mem back there.(!x) then
          Hashtbl.find back there.(!x) else -1
    done;
    while !start < alen & used.(!start) do start := !start + 1 done;
    if !start = alen then start := -1
  done;
  permute mixed (invperm perm);;

(* works only for limited length arrays 2^10 *)
let asample k a = atake k (Random.shuffle (Array.enum a));;


(* stochastic universal sampling *)
let sus (randval:float) (n:int) (pop:popt) =
  let sw = List.fold_left (+.)  0.0 (List.map fst pop) in
  let inc = sw/.(float n) in
  let pairs = List.rev (List.sort ~cmp:Pervasives.compare pop) in
  let rec next p t i = match p with
      [] -> failwith "sus panic"
    | ((xw,x)::xs) ->
      if t -. xw < 0.0 then (p, t +. inc, i, (xw,x))
      else next xs (t -. xw) (i + 1)
  in
  let rec iter p t i n acc = if n < 1 then acc else
      let (p', t', i', x') = next p t i in iter p' t' i' (n - 1) (x'::acc) 
  in List.rev (iter pairs (randval/.(float n) *. sw) 0 n []);;



(* generates a function that performs xover *)
let makexo (fitness:fitt) =
  (fun (af,a) (bf,b) ->
    List.map (fun v -> (fitness v, v)) [(mix a b); (mix a b)]);;

(* generates a function that mutatates *)
let makemut (fitness:fitt) (ri:unit -> indt) =
  (fun (vf,v) -> let w = mix v (snd (ri ())) in (fitness w, w));;

(* quick and blah shuffle of lists... *)
let blahshuffle l =
  let cmp _ _ = (Random.int 2) * 2 - 1 in List.sort ~cmp:cmp l;;

(* sample a quarter of the population to be parents and mate them,
   return list of offspring *)
let crossover (xo:indt->indt->indt list) (pop:popt) =
  let plen = List.length pop in
  let parents = sus (Random.float 1.) (plen/4) pop in
  let (moms, dads) = List.split_at (plen/8) (blahshuffle parents) in
  List.concat (List.map2 xo moms dads);;

(* sample psize/div individuals and mutate,
   returning results in a list *)
let mutate (mut:indt -> indt) (div:int) (pop:popt) =
  let plen = List.length pop in
  List.map mut (sus (Random.float 1.) (plen/div) pop);;

(* loop a function n times  returning results in a list *)
let gloop g init n =
  let rec helper g n y acc =
    if n > 0 then
      let x = g y in helper g (n - 1) x (x::acc)
    else acc in
  List.rev (helper g n init []);;

(* remove duplicates and replenish with random elements *)
let diversify ri (pop:popt) =
  let n = List.length pop and
      newpop = List.sort_unique Pervasives.compare pop in
  if List.length newpop = n then newpop else
    newpop @ (gloop (fun _ -> ri()) (0.,[|1|]) (n - (List.length newpop)));;

(* create an initial random population, using ri to generate each
   individual. Return a list of n individuals. *)
let initpop ri n =
  let il = List.of_enum (1--n) in List.map (fun _ -> ri ()) il;;


(* do one generation: generational and elitist *)
let generation (xo:indt->indt->indt list) (mut:indt->indt) (diversify:popt->popt) (div:int) (pop:popt) =
  let ox = crossover xo pop and
      om = mutate mut div pop in
  let cand = List.rev (List.sort (diversify (List.concat [om ; ox ; pop]))) in
  (List.hd cand) :: (sus (Random.float 1.) ((List.length pop) - 1) cand) ;;

(* run the algorithm, using a function gen : pop -> newpop.
   Waits wait generations without improvement before halting. *)
let run ?l:((log, strm) = (false, stderr)) wait gen pop = 
  let rec iter delay pop (bestw,bestx) =
    if delay < 1 then (bestw,bestx) else
      let ((xw,x)::rest) = gen pop in
      if log then
        IO.write_line strm (String.of_float xw);
      if xw > bestw then iter wait ((xw,x)::rest) (xw,x)
      else iter (delay - 1) ((xw,x)::rest) (bestw,bestx)
  in
  iter wait pop (List.max pop);;

(* the entry point to the algorithm for the assignment problem *)
let runaga ?l:((log, strm) = (false, stderr)) wait psize ri fit =
  let xo = makexo fit and
      mut = makemut fit ri in
  (* let gen = generation xo mut (diversify ri) 10 in *)
  let gen = generation xo mut (fun x -> x) 10 in
  run ~l:(log,strm) wait gen (initpop ri psize);;
