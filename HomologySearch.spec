/*
A KBase module: HomologySearch
*/

module HomologySearch {
  /*
    KBase genome object id.
  */

  typedef structure {
    string workspace;
    string sequence;
    string database;
    string search_type;
    list<string> genome_ids;
    string program;
    int max_hits;
    string evalue_cutoff;
  } HomologySearchInputParams;

  /*
    output format from
    https://github.com/kbase/genome_util/blob/master/KBaseGenomeUtil.spec
  */
  typedef structure {
    string Parameters_expect;
    string Parameters_filter;
    string Parameters_gap-extend;
    string Parameters_gap-open;
    string Parameters_matrix;
    string Parameters_sc-match;
    string Parameters_sc-mismatch;
  } Parameters;

  typedef structure {
    Parameters Parameters;
  } BlastOutput_param;

  typedef structure {
    string Hsp_align-len;
    string Hsp_bit-score;
    string Hsp_evalue;
    string Hsp_hit-frame;
    string Hsp_hit-from;
    string Hsp_hit-to;
    string Hsp_hseq;
    string Hsp_identity;
    string Hsp_midline;
    string Hsp_num;
    string Hsp_positive;
    string Hsp_qseq;
    string Hsp_query-frame;
    string Hsp_query-from;
    string Hsp_query-to;
    string Hsp_score;
  } Hsp_details;

  typedef list <Hsp_details> Hsp;

  typedef structure {
    Hsp Hsp;
  } Hit_hsps;

  typedef structure {
    string Hit_accession;
    string Hit_def;
    string Hit_id;
    string Hit_len;
    string Hit_num;
    Hit_hsps Hit_hsps;
  } hit_details;

  typedef list <hit_details> Hit;

  typedef structure {
    Hit Hit;
  } Iteration_hits;

  typedef structure {
    Iteration_hits Iteration_hits;
    string Iteration_iter-num;
    string Iteration_query-ID;
    string Iteration_query-def;
    string Iteration_query-len;
  } Iteration_details;

  typedef list <Iteration_details> Iteration;

  typedef structure {
    Iteration Iteration;
  } BlastOutput_iterations;

  typedef structure {
    string BlastOutput_db;
    string BlastOutput_program;
    string BlastOutput_query-ID;
    string BlastOutput_query-def;
    string BlastOutput_query-len;
    string BlastOutput_reference;
    string BlastOutput_version;
    BlastOutput_param BlastOutput_param;
    BlastOutput_iterations BlastOutput_iterations;
  } BlastOutput;

  typedef structure {
    string json_output;
  } BlastJSONOutput;

  /*
  methods
  */
  funcdef blast_fasta(HomologySearchInputParams params)
    /* returns (BlastOutput) */
    returns (BlastJSONOutput output)
    authentication required;
};
