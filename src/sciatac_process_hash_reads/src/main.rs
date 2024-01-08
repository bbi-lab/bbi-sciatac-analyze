



/*
** Turn off warnings about unused parentheses.
*/
#![allow(unused_parens)]


use argparse::ArgumentParser;
use std::fs;
use std::io;
use std::io::BufRead;
use std::io::Write;


/*
** Input file columns
**   col #   description
**   -----   -----------
**    1       sciPlexATAC    (string)
**    2       cell_id        (string)
**    3       umi            (string)
**    4       hash name      (string)
**    5       hash axis?     (integer)
*/
fn read_hash_read_file(input_filename: &str) -> io::Result<Vec<String>> {
  let fp = fs::File::open(input_filename)?;
  let bfr = io::BufReader::new(fp);
  bfr.lines().collect::<io::Result<Vec<String>>>()
}


/*
** Shift UMI token to the end of the strings.
** Note: this enables simple filtering after
** sorting.
*/
fn reform_hash_read_lines_01(vec_in: Vec<String>) -> io::Result<Vec<String>> {
  let mut vec_out: Vec<String> = vec![];

  for line in vec_in {
    let toks: Vec<&str> = line.split('\t').collect();
    vec_out.push([toks[0], toks[1], toks[3], toks[4], toks[2]].join("\t"));
  }
  Ok(vec_out)
}


/*
** Shift first token to the end of the strings.
** Note: this enables simple filtering after
** sorting.
*/
fn reform_hash_read_lines_02(vec_in: Vec<String>) -> io::Result<Vec<String>> {
  let mut vec_out: Vec<String> = vec![];

  for line in vec_in {
    let toks: Vec<&str> = line.split('\t').collect();
    vec_out.push([toks[1], toks[2], toks[3], toks[4], toks[0]].join("\t"));
  }
  Ok(vec_out)
}


fn main() -> Result<(), Box<dyn std::error::Error>> {

  let mut input_filename = String::new();
  let mut sample_name    = String::new();

  /*
  ** Note: try 'clap' command line argument package.
  */
  {
    let mut parser = ArgumentParser::new();
    parser.set_description("Process hash reads from table made by sciatac_find_hash_reads.");

    parser.refer(&mut input_filename).required().add_option(&["-i", "--input_filename"], argparse::Store, "Name of hash record file.");
    parser.refer(&mut sample_name).required().add_option(&["-s", "--sample_name"], argparse::Store, "Name of sample.");
    parser.add_option(&["-V", "--version"], argparse::Print(env!("CARGO_PKG_VERSION").to_string()), "Show version.");
    parser.parse_args_or_exit();
  }

  if(input_filename.is_empty() ||
     sample_name.is_empty()) {
    eprintln!("sciatac_process_hash_reads: error: missing command line argument\nUsage: sciatac_process_hash_reads -i <hash_record_filename> -s <sample_name>");
    std::process::exit(1);
  }

  /*
  ** Make hashTable.out, which has the distinct counts
  ** by cell + hash id (and columns 1 and 5 of the
  ** original file).
  */
  let hash_table_out_filename = format!("{}-{}", sample_name, "hashTable.out");
  let out_file = fs::File::create(hash_table_out_filename).expect("error: bad status: fs::File::Create");
  let mut buf_out_file = io::BufWriter::with_capacity(1024*1024, out_file);

  let mut vec_lines = read_hash_read_file(&input_filename).expect("error: read_hash_record_table");
  vec_lines = reform_hash_read_lines_01(vec_lines).expect("error: reform_hash_read_lines_01");
  vec_lines.sort();

  let mut line_prv: String = String::new();
  let mut toks_prv: Vec<&str> = vec!["", "", "", "", ""];
  let mut toks_cur: Vec<&str>;

  let mut n1: u64 = 0;

  for line_cur in vec_lines {

    /*
    ** Drop identical lines, which are considered to be PCR duplicates.
    */
    if(line_cur == line_prv) {
      continue;
    }

    toks_cur = line_cur.split('\t').collect();

    /*
    ** Report cell + hash combinations (and tokens 0 and 3).
    */
    if(n1 > 0 &&
       (toks_cur[0] != toks_prv[0] ||
        toks_cur[1] != toks_prv[1] ||
        toks_cur[2] != toks_prv[2] ||
        toks_cur[3] != toks_prv[3])) {
      writeln!(buf_out_file, "{}\t{}\t{}\t{}\t{}", toks_prv[0], toks_prv[1], toks_prv[2], toks_prv[3], n1)?;
      n1 = 0;
    }

    /*
    ** The UMI sequence has changed. (All other tokens
    ** must be the same as in the previous line and
    **  identical lines are skipped.)
    */
    n1 += 1;

    /*
    ** Prepare for the next loop iteration.
    */
    line_prv = line_cur.to_string();
    toks_prv = line_prv.split('\t').collect();
  }

  /*
  ** Report the last entry.
  */
  if(n1 > 0) {
    writeln!(buf_out_file, "{}\t{}\t{}\t{}\t{}", toks_prv[0], toks_prv[1], toks_prv[2], toks_prv[3], n1)?;
  }

  buf_out_file.flush().expect("error: bad status: BufWriter.flush()");


  /*
  ********************************************************************************* 
  */


  /*
  ** Make hashReads.per.cell and hashUMIs.per.cell.
  */
  let hash_reads_out_filename = format!("{}-{}", sample_name, "hashReads.per.cell");
  let out_file_1 = fs::File::create(hash_reads_out_filename).expect("error: bad status: fs::File::Create");
  let mut buf_out_file_1 = io::BufWriter::with_capacity(1024*1024, out_file_1);

  let hash_umis_out_filename = format!("{}-{}", sample_name, "hashUMIs.per.cell");
  let out_file_2 = fs::File::create(hash_umis_out_filename).expect("error: bad status: fs::File::Create");
  let mut buf_out_file_2 = io::BufWriter::with_capacity(1024*1024, out_file_2);

  let mut vec_lines = read_hash_read_file(&input_filename).expect("error: read_hash_record_table");
  vec_lines = reform_hash_read_lines_02(vec_lines).expect("error: reform_hash_read_lines_02");
  vec_lines.sort();

  let mut line_prv: String = String::new();
  let mut toks_prv: Vec<&str> = vec!["", "", "", "", ""];
  let mut toks_cur: Vec<&str>;

  let mut n1: u64 = 0;
  let mut n2: u64 = 0;

  for line_cur in vec_lines {

    toks_cur = line_cur.split('\t').collect();

    /* 
    ** Report reads per cell.
    */
    if(n1 > 0 &&
       (toks_cur[0] != toks_prv[0])) {
      writeln!(buf_out_file_1, "{}\t{}\t{}", sample_name, toks_prv[0], n1)?;
      n1 = 0;
    }

    /*
    ** Count reads per cell.
    */
    n1 += 1;

    /*
    ** Drop identical lines, which are considered to be PCR duplicates.
    */
    if(line_cur == line_prv) {
      continue;
    }

    /*
    ** Report UMIs per cell.
    */
    if(n2 > 0 &&
       (toks_cur[0] != toks_prv[0])) {
      writeln!(buf_out_file_2, "{}\t{}\t{}", sample_name, toks_prv[0], n2)?;
      n2 = 0;
    }

    /*
    ** Count UMIs per cell.
    */
    n2 += 1;

    /*
    ** Prepare for the next loop iteration.
    */
    line_prv = line_cur.to_string();
    toks_prv = line_prv.split('\t').collect();
  }

  /*
  ** Report the last entry.
  */
  if(n1 > 0) {
    writeln!(buf_out_file_1, "{}\t{}\t{}", sample_name, toks_prv[0], n1)?;
  }
  buf_out_file_1.flush().expect("error: bad status: BufWriter.flush()");


  if(n2 > 0) {
    writeln!(buf_out_file_2, "{}\t{}\t{}", sample_name, toks_prv[0], n2)?;
  }
  buf_out_file_2.flush().expect("error: bad status: BufWriter.flush()");


 Ok(())
}
