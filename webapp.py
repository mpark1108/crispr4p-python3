#!/usr/bin/env python3
# Standalone local HTTP web server for CRISPR4P (Python 3 compatible)

import os
import urllib.parse
import sys
from http.server import HTTPServer, BaseHTTPRequestHandler

from crispr4p.service import Crispr4pService, OligoLengthError
from crispr4p.web_views import (
    build_spedit_candidate_data,
    render_design_result,
    render_execution_error,
    render_missing_query_error,
    render_oligo_length_error,
    render_oligo_result,
)


PORT = 8080


def create_service():
    """Create an isolated application service for one web operation."""
    return Crispr4pService.from_project_data(
        precomputed_folder="precomputed",
    )


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
        try:
            result = create_service().analyze_oligo(
                oligo_seq,
                n_mismatch=mismatches,
            )
        except OligoLengthError as error:
            return render_oligo_length_error(error.sequence_length)

        return render_oligo_result(result)

    def run_design_model(self, name, chromosome, coor_lower, coor_upper):
        service = create_service()
        if name is not None:
            result = service.design_gene(name, n_mismatch=0)
        else:
            result = service.design_region(
                chromosome,
                coor_lower,
                coor_upper,
                n_mismatch=0,
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
