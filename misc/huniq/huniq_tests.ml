open Batteries
open Printf
open Should
open OUnit

Printexc.record_backtrace true

let open_huniq args =
    let fn_exe = Filename.concat (Filename.dirname Sys.argv.(0)) "huniq.native"
    Unix.open_process (fn_exe ^ " " ^ (String.concat " " args))

let simple_test args input_lines expected_lines =
    let huniq_in, huniq_out = open_huniq args
    foreach input_lines
        fun line ->
            ignore (IO.really_output huniq_out line 0 (String.length line))
            IO.write huniq_out '\n'
    close_out huniq_out
    (huniq_in |> IO.lines_of |> List.of_enum) $hould # equal expected_lines

let trivial1 () = simple_test [] (List.enum ["foo"]) ["foo"]
let trivial2 () = simple_test [] (List.enum ["foo"; "foo"]) ["foo"]
let trivial3 () = simple_test [] (List.enum ["foo"; "bar"]) ["foo"; "bar"]
let trivial4 () = simple_test [] (List.enum ["foo"; "bar"; "foo"]) ["foo"; "bar"]
let trivial5 () = simple_test [] (List.enum ["foo"; "bar"; "foo"; "bar"]) ["foo"; "bar"]
let trivial6 () = simple_test [] (List.enum ["foo"; "bar"; "bar"; "foo"]) ["foo"; "bar"]

let count1 () = simple_test ["-c"] (List.enum ["foo"]) ["1 foo"]
let count2 () = simple_test ["-c"] (List.enum ["foo"; "foo"]) ["2 foo"]
let count3 () = simple_test ["-c"] (List.enum ["foo"; "bar"]) ["1 foo"; "1 bar"]
let count4 () = simple_test ["-c"] (List.enum ["foo"; "bar"; "foo"]) ["2 foo"; "1 bar"]
let count5 () = simple_test ["-c"] (List.enum ["foo"; "bar"; "foo"; "bar"]) ["2 foo"; "2 bar"]
let count6 () = simple_test ["-c"] (List.enum ["foo"; "bar"; "bar"; "foo"]) ["2 foo"; "2 bar"]

let randomized () =
    for i = 1 to 100 do
        let i1 = Random.int 1000
        let i2 = Random.int 1000
        let i3 = Random.int 1000
        let ar = Array.concat [Array.make i1 "foo"; Array.make i2 "bar"; Array.make i3 "ba z"]
        let ar = Random.shuffle (Array.enum ar)
        let ar = Array.concat [[|"foo";"bar";"ba z"|]; ar]

        simple_test [""] (Array.enum ar) ["foo";"bar";"ba z"]
        simple_test ["-c"] (Array.enum ar) [sprintf "%d foo" (i1+1);
                                            sprintf "%d bar" (i2+1);
                                            sprintf "%d ba z" (i3+1)]


let all_tests = "should tests" >::: [
        "trivial1" >:: trivial1;
        "trivial2" >:: trivial2;
        "trivial3" >:: trivial3;
        "trivial4" >:: trivial4;
        "trivial5" >:: trivial5;
        "trivial6" >:: trivial6;
        "count1" >:: count1;
        "count2" >:: count2;
        "count3" >:: count3;
        "count4" >:: count4;
        "count5" >:: count5;
        "count6" >:: count6;
        "randomized" >:: randomized
    ]

Random.init 12345
run_test_tt_main all_tests
