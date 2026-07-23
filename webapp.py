#!/usr/bin/env python3
# Standalone local HTTP web server for CRISPR4P (Python 3 compatible)

import html
import os
import json
import urllib.parse
import sys
from http.server import HTTPServer, BaseHTTPRequestHandler

import crispr4p.crispr4p as crp
from crispr4p.spedit import (
    has_internal_bsai_site,
    make_spedit_oligos,
)


PORT = 8080


def build_spedit_candidate_data(table_pos_grna) -> list[dict]:
    """
    Build SpEDIT output aligned one-to-one with tablePos_grna.

    The returned list uses the same candidate ordering and index as the
    existing browser-side guide and primer selector.
    """
    candidates = []

    for row in table_pos_grna:
        try:
            # Existing CRISPR4P result layout:
            # row[1][0] is the complete 20-nt guide.
            guide = row[1][0]
            forward, reverse = make_spedit_oligos(guide)

            candidates.append(
                {
                    "guide": guide,
                    "forward": forward,
                    "reverse": reverse,
                    "has_internal_bsai": has_internal_bsai_site(guide),
                    "error": None,
                }
            )

        except (IndexError, TypeError, ValueError) as error:
            candidates.append(
                {
                    "guide": "",
                    "forward": "",
                    "reverse": "",
                    "has_internal_bsai": False,
                    "error": str(error),
                }
            )

    return candidates


class CRISPR4PHandler(BaseHTTPRequestHandler):
    def do_GET(self):
        parsed_path = urllib.parse.urlparse(self.path)
        path = parsed_path.path

        if path == "/" or path == "/index.html" or path == "/webapp.py":
            self.serve_form()
        elif path.startswith("/css/"):
            self.serve_css(path)
        else:
            self.send_error(404, "File not found")

    def do_POST(self):
        parsed_path = urllib.parse.urlparse(self.path)
        path = parsed_path.path

        if path == "/webapp.py" or path == "/":
            self.process_post()
        else:
            self.send_error(404, "File not found")

    def serve_form(self, result_content=""):
        src_path = os.path.dirname(__file__)
        src_path = "." if src_path == "" else src_path

        try:
            template_path = os.path.join(src_path, 'template/bahler_template.html')
            with open(template_path, 'r', encoding='utf-8') as fh:
                template_file = fh.read()
            
            # The template has "%s" at line 61 which is replaced by result_content
            rendered = template_file % result_content
            
            self.send_response(200)
            self.send_header("Content-type", "text/html; charset=utf-8")
            self.end_headers()
            self.wfile.write(rendered.encode('utf-8'))
        except IOError as err:
            self.send_error(500, f"Error reading templates: {err}")

    def serve_css(self, path):
        src_path = os.path.dirname(__file__)
        src_path = "." if src_path == "" else src_path
        filename = os.path.basename(path)
        css_file = os.path.join(src_path, "css", filename)
        if os.path.exists(css_file):
            try:
                with open(css_file, 'r', encoding='utf-8') as f:
                    content = f.read()
                self.send_response(200)
                self.send_header("Content-type", "text/css")
                self.end_headers()
                self.wfile.write(content.encode('utf-8'))
            except IOError:
                self.send_error(500, "Error reading CSS file")
        else:
            self.send_error(404, "CSS File not found")

    def process_post(self):
        content_length = int(self.headers.get('Content-Length', 0))
        post_data = self.rfile.read(content_length).decode('utf-8')
        params = urllib.parse.parse_qs(post_data)

        # Extract values
        name = params.get('name', [None])[0]
        chromosome = params.get('chromosome', [None])[0]
        coor_lower = params.get('coor_lower', [None])[0]
        coor_upper = params.get('coor_upper', [None])[0]
        oligo_sequence = params.get('oligo_sequence', [None])[0]
        oligo_mismatch_str = params.get('oligo_mismatch', ['0'])[0]

        # Clean values
        name = name.strip() if name else None
        chromosome = chromosome.strip() if chromosome else None
        coor_lower = coor_lower.strip() if coor_lower else None
        coor_upper = coor_upper.strip() if coor_upper else None
        oligo_sequence = oligo_sequence.strip().upper() if oligo_sequence else None

        result_html = ""

        try:
            if oligo_sequence:
                try:
                    mismatches = int(oligo_mismatch_str)
                except ValueError:
                    mismatches = 0
                result_html = self.run_oligo_model(oligo_sequence, mismatches)
            elif name or (chromosome and coor_lower and coor_upper):
                result_html = self.run_design_model(name, chromosome, coor_lower, coor_upper)
            else:
                result_html = '<font color="red"><h3>Error: Please fill either Name, Coordinates, or Oligo Sequence</h3></font>'
        except Exception as e:
            result_html = f'<font color="red"><h3>ERROR during execution: {str(e)}</h3></font>'

        self.serve_form(result_html)

    def run_oligo_model(self, oligo_seq, mismatches):
        if len(oligo_seq) == 20:
            seed = oligo_seq
            pam = "NGG"
        elif len(oligo_seq) == 23:
            seed = oligo_seq[:20]
            pam = oligo_seq[20:]
        else:
            return f'<font color="red"><h3>Error: Oligo sequence must be 20 bp (seed only) or 23 bp (seed + PAM). Current length: {len(oligo_seq)}</h3></font>'

        spedit_forward, spedit_reverse = make_spedit_oligos(seed)

        if has_internal_bsai_site(seed):
            spedit_warning = (
                '<strong style="color: #b00020;">'
                "Warning: this guide contains an internal BsaI recognition site."
                "</strong>"
            )
        else:
            spedit_warning = "No internal BsaI site detected."

        datapath = "data/"
        FASTA = datapath + 'Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa'
        COORDINATES = datapath + 'COORDINATES.txt'
        SYNONIMS = datapath + 'SYNONIMS.txt'

        pd = crp.PrimerDesign(FASTA, COORDINATES, SYNONIMS, precomputed_folder='precomputed')
        pd.getNGGsFromGenome()

        query_ngg = crp.NGG(chro='query', pos=0, strand=1, seed=seed, pam=pam)
        _, tableDict = pd._single_table_worker(query_ngg, mismatches)

        match_8 = len(tableDict.get(8, []))
        match_10 = len(tableDict.get(10, []))
        match_12 = len(tableDict.get(12, []))
        match_14 = len(tableDict.get(14, []))
        match_16 = len(tableDict.get(16, []))
        match_18 = len(tableDict.get(18, []))
        match_20 = len(tableDict.get(20, []))

        matches_20 = tableDict.get(20, [])
        details_html = ""
        if matches_20:
            details_html += '<table border="1" cellpadding="5" cellspacing="0" style="border-collapse: collapse; width: 100%; font-family: monospace; font-size: 12px;">'
            details_html += '<tr style="background-color: #D1F0A6;"><th>#</th><th>Chromosome</th><th>Position</th><th>Strand</th><th>Genomic Target Sequence (Seed)</th><th>PAM</th></tr>'
            for idx, match in enumerate(matches_20):
                strand_str = "+" if match.strand == 1 else "-"
                details_html += f'<tr><td>{idx+1}</td><td>{match.chromosome}</td><td>{match.pos}</td><td>{strand_str}</td><td>{match.seed}</td><td>{match.pam}</td></tr>'
            details_html += '</table>'
        else:
            details_html = "<p>No full 20bp target/off-target matches found in the genome.</p>"

        block = f"""
        <div id="search_content">
          <div id="search_summary">
              <h4>Oligo Search Results:</h4>
              <b>Oligo Sequence (Query)</b>: {oligo_seq}<br>
              <b>Seed Segment (20bp)</b>: {seed}<br>
              <b>Mismatches Allowed</b>: {mismatches}<br>
              <hr>

              <h4>SpEDIT/pLSB BsaI Golden Gate oligos</h4>

              <b>Forward oligo, 52 nt, 5&#8242;&rarr;3&#8242;</b>:
              <code>{spedit_forward}</code><br>

              <b>Reverse oligo, 52 nt, 5&#8242;&rarr;3&#8242;</b>:
              <code>{spedit_reverse}</code><br>

              <b>Internal BsaI site check</b>: {spedit_warning}<br>
          </div>

          <h3 class="toggle_header">Genome Match Summary</h3>
          <div style="padding: 10px;">
            <table border="1" cellpadding="5" cellspacing="0" style="border-collapse: collapse; width: 50%;">
              <thead>
                <tr style="background-color: #D1F0A6;">
                  <th>Seed Match Length</th>
                  <th>Matching Sites (adjacent to NGG/NAG PAM)</th>
                </tr>
              </thead>
              <tbody>
                <tr><td>8 bp</td><td>{match_8}</td></tr>
                <tr><td>10 bp</td><td>{match_10}</td></tr>
                <tr><td>12 bp</td><td>{match_12}</td></tr>
                <tr><td>14 bp</td><td>{match_14}</td></tr>
                <tr><td>16 bp</td><td>{match_16}</td></tr>
                <tr><td>18 bp</td><td>{match_18}</td></tr>
                <tr><td>20 bp</td><td>{match_20}</td></tr>
              </tbody>
            </table>
          </div>

          <h3 class="toggle_header">Details of Full 20bp Matches</h3>
          <div style="padding: 10px;">
            {details_html}
          </div>
        </div>
        """
        return block

    def run_design_model(self, name, chromosome, coor_lower, coor_upper):
        datapath = "data/"
        FASTA = datapath + 'Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa'
        COORDINATES = datapath + 'COORDINATES.txt'
        SYNONIMS = datapath + 'SYNONIMS.txt'

        pd = crp.PrimerDesign(FASTA, COORDINATES, SYNONIMS, precomputed_folder='precomputed')
        tablePos_grna, hr_dna, primercheck, name, cr, start, end = pd.runWeb(name, chromosome, coor_lower, coor_upper, nMismatch=0)

        pm = primercheck[0] if primercheck else {}
        
        def get_tm_str(val):
            try:
                return "%d &deg;C" % int(round(float(val)))
            except (ValueError, TypeError):
                return "- &deg;C"

        result_dict = {'name': name or '-',
                       'chromosome': cr,
                       'start': start,
                       'end': end,
                       'hrfw': hr_dna[0],
                       'hrrv': hr_dna[1],
                       'deleted_dna': hr_dna[2],
                       'primer_left': pm.get('PRIMER_LEFT_0_SEQUENCE', '-'),
                       'left_tm': get_tm_str(pm.get('PRIMER_LEFT_0_TM', 0)),
                       'primer_right': pm.get('PRIMER_RIGHT_0_SEQUENCE', '-'),
                       'right_tm': get_tm_str(pm.get('PRIMER_RIGHT_0_TM', 0)),
                       'deleted_dna_size': str(pm.get('PRIMER_PAIR_0_PRODUCT_SIZE', '-')) + " (bp)",
                       'negative_result_size': str(pm.get('negative_result', '-')) + " (bp)"}

        result_dict['json_table'] = json.dumps(tablePos_grna)
        result_dict["spedit_json"] = json.dumps(
            build_spedit_candidate_data(tablePos_grna)
        )

        src_path = os.path.dirname(__file__) if os.path.dirname(__file__) else '.'
        with open(os.path.join(src_path, 'template/container_table.html'), 'r', encoding='utf-8') as fh:
            template_file = fh.read()

        return template_file % (result_dict)


def main():
    server_address = ('', PORT)
    httpd = HTTPServer(server_address, CRISPR4PHandler)
    print(f"Starting local server on http://localhost:{PORT} ...")
    print("Press Ctrl+C to stop.")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down server.")
        sys.exit(0)


if __name__ == "__main__":
    main()
