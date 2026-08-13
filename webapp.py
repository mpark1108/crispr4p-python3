#!/usr/bin/env python3
# Standalone local HTTP web server for CRISPR4P (Python 3 compatible)

import os
import urllib.parse
import sys
from http.server import HTTPServer, BaseHTTPRequestHandler

import crispr4p.crispr4p as crp
from crispr4p.models import DesignResult
from crispr4p.spedit import (
    has_internal_bsai_site,
    make_spedit_oligos,
)
from crispr4p.web_views import (
    OligoMatchView,
    build_spedit_candidate_data,
    render_design_result,
    render_execution_error,
    render_missing_query_error,
    render_oligo_length_error,
    render_oligo_result,
)


PORT = 8080


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
                result_html = render_missing_query_error()
        except Exception as e:
            result_html = render_execution_error(e)

        self.serve_form(result_html)

    def run_oligo_model(self, oligo_seq, mismatches):
        if len(oligo_seq) == 20:
            seed = oligo_seq
            pam = "NGG"
        elif len(oligo_seq) == 23:
            seed = oligo_seq[:20]
            pam = oligo_seq[20:]
        else:
            return render_oligo_length_error(len(oligo_seq))

        spedit_forward, spedit_reverse = make_spedit_oligos(seed)

        datapath = "data/"
        FASTA = datapath + 'Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa'
        COORDINATES = datapath + 'COORDINATES.txt'
        SYNONIMS = datapath + 'SYNONIMS.txt'

        pd = crp.PrimerDesign(FASTA, COORDINATES, SYNONIMS, precomputed_folder='precomputed')
        pd.getNGGsFromGenome()

        query_ngg = crp.NGG(chro='query', pos=0, strand=1, seed=seed, pam=pam)
        _, tableDict = pd._single_table_worker(query_ngg, mismatches)

        match_counts = {
            length: len(tableDict.get(length, []))
            for length in (8, 10, 12, 14, 16, 18, 20)
        }
        full_matches = []
        for match in tableDict.get(20, []):
            pam_coordinates, cut_coordinates = pd.getOligoHitCoordinates(match)
            full_matches.append(
                OligoMatchView(
                    chromosome=match.chromosome,
                    pam_coordinates=pam_coordinates,
                    cut_coordinates=cut_coordinates,
                    strand=match.strand,
                    seed=match.seed,
                    pam=match.pam,
                )
            )

        return render_oligo_result(
            oligo_sequence=oligo_seq,
            seed=seed,
            mismatches=mismatches,
            spedit_forward=spedit_forward,
            spedit_reverse=spedit_reverse,
            has_internal_bsai=has_internal_bsai_site(seed),
            match_counts=match_counts,
            full_matches=full_matches,
        )

    def run_design_model(self, name, chromosome, coor_lower, coor_upper):
        datapath = "data/"
        FASTA = datapath + 'Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa'
        COORDINATES = datapath + 'COORDINATES.txt'
        SYNONIMS = datapath + 'SYNONIMS.txt'

        pd = crp.PrimerDesign(FASTA, COORDINATES, SYNONIMS, precomputed_folder='precomputed')
        result = DesignResult.from_legacy(
            pd.runWeb(
                name,
                chromosome,
                coor_lower,
                coor_upper,
                nMismatch=0,
            )
        )

        src_path = os.path.dirname(__file__) if os.path.dirname(__file__) else '.'
        with open(os.path.join(src_path, 'template/container_table.html'), 'r', encoding='utf-8') as fh:
            template_file = fh.read()

        return render_design_result(result, template_file)


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
