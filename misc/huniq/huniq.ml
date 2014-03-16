open Batteries
open Printf
open Extlib.OptParse

let usage = "%prog [OPTION]... [INPUT [OUTPUT]]

Filter repeated lines from INPUT (or standard input),
writing to OUTPUT (or standard output). Lines need not
be in sorted order, and memory usage is proportional
to the total unique input lines."

let opt_parser =
    OptParser.make
        ~usage
        ()

let opt_count = StdOpt.store_true ()
OptParser.add ~help:"prefix lines by the number of occurrences" ~short_name:'c' ~long_name:"count" opt_parser opt_count

let pos_args = Array.of_list (OptParser.parse_argv opt_parser)

if Array.length pos_args > 2 then
    OptParser.usage opt_parser ()
    exit 1

let opt_count = Opt.get opt_count

let infile =
    if Array.length pos_args = 0 || pos_args.(0) = "-" then stdin
    else File.open_in ~mode:[`text] pos_args.(0)

let outfile =
    if Array.length pos_args < 2 || pos_args.(1) = "-" then stdout
    else File.open_out ~mode:[`create;`text;`trunc] pos_args.(1)

let main _ =
    let tbl = Hashtbl.create 64
    let ord = Queue.create ()
    foreach (IO.lines_of infile)
        fun line ->
            let occ = Hashtbl.find_default tbl line 0
            if occ = 0 then
                if opt_count then Queue.add line ord
                else IO.write_line outfile line
            Hashtbl.replace tbl line (occ+1)
    if opt_count then
        foreach (Queue.enum ord)
            fun line -> fprintf outfile "%d %s\n" (Hashtbl.find tbl line) line

with_dispose ~dispose:flush main outfile
