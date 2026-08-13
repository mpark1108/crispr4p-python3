"""Pure HTML rendering helpers for CRISPR4P web results."""

import json

from .models import DesignResult, OligoAnalysisResult, OligoMatch
from .spedit import has_internal_bsai_site, make_spedit_oligos


# Retain the step-4 import name while the service model becomes canonical.
OligoMatchView = OligoMatch


def build_spedit_candidate_data(table_pos_grna) -> list[dict]:
    """Build SpEDIT output aligned one-to-one with the legacy guide table."""
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


def render_missing_query_error():
    return (
        '<font color="red"><h3>Error: Please fill either Name, '
        'Coordinates, or Oligo Sequence</h3></font>'
    )


def render_execution_error(error):
    return f'<font color="red"><h3>ERROR during execution: {str(error)}</h3></font>'


def render_oligo_length_error(sequence_length):
    return (
        '<font color="red"><h3>Error: Oligo sequence must be 20 bp '
        f'(seed only) or 23 bp (seed + PAM). Current length: {sequence_length}'
        '</h3></font>'
    )


def render_oligo_result(result: OligoAnalysisResult):
    """Render an oligo analysis without accessing HTTP or genome state."""
    if result.has_internal_bsai:
        spedit_warning = (
            '<strong style="color: #b00020;">'
            "Warning: this guide contains an internal BsaI recognition site."
            "</strong>"
        )
    else:
        spedit_warning = "No internal BsaI site detected."

    details_html = ""
    if result.full_matches:
        details_html += '<table border="1" cellpadding="5" cellspacing="0" style="border-collapse: collapse; width: 100%; font-family: monospace; font-size: 12px;">'
        details_html += '<tr style="background-color: #D1F0A6;"><th>#</th><th>Chromosome</th><th>PAM coordinates (1-based, inclusive)</th><th>Cas9 cut</th><th>Strand</th><th>Genomic Target Sequence (Seed)</th><th>PAM</th></tr>'
        for index, match in enumerate(result.full_matches):
            strand = "+" if match.strand == 1 else "-"
            details_html += (
                f'<tr><td>{index+1}</td><td>{match.chromosome}</td>'
                f'<td>{match.pam_coordinates[0]} - {match.pam_coordinates[1]}</td>'
                f'<td>{match.cut_coordinates[0]} | {match.cut_coordinates[1]}</td>'
                f'<td>{strand}</td><td>{match.seed}</td>'
                f'<td>{match.pam}</td></tr>'
            )
        details_html += '</table>'
    else:
        details_html = "<p>No full 20bp target/off-target matches found in the genome.</p>"

    block = f"""
        <div id="search_content">
          <div id="search_summary">
              <h4>Oligo Search Results:</h4>
              <b>Oligo Sequence (Query)</b>: {result.oligo_sequence}<br>
              <b>Seed Segment (20bp)</b>: {result.seed}<br>
              <b>Mismatches Allowed</b>: {result.n_mismatch}<br>
              <hr>

              <h4>SpEDIT/pLSB BsaI Golden Gate oligos</h4>

              <b>Forward oligo, 52 nt, 5&#8242;&rarr;3&#8242;</b>:
              <code>{result.spedit_forward}</code><br>

              <b>Reverse oligo, 52 nt, 5&#8242;&rarr;3&#8242;</b>:
              <code>{result.spedit_reverse}</code><br>

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
                <tr><td>8 bp</td><td>{result.match_counts.get(8, 0)}</td></tr>
                <tr><td>10 bp</td><td>{result.match_counts.get(10, 0)}</td></tr>
                <tr><td>12 bp</td><td>{result.match_counts.get(12, 0)}</td></tr>
                <tr><td>14 bp</td><td>{result.match_counts.get(14, 0)}</td></tr>
                <tr><td>16 bp</td><td>{result.match_counts.get(16, 0)}</td></tr>
                <tr><td>18 bp</td><td>{result.match_counts.get(18, 0)}</td></tr>
                <tr><td>20 bp</td><td>{result.match_counts.get(20, 0)}</td></tr>
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


def render_design_result(result: DesignResult, template_text):
    """Render a structured design result using the existing HTML template."""
    primer = result.checking_primers[0] if result.checking_primers else {}

    def get_tm_str(value):
        try:
            return "%d &deg;C" % int(round(float(value)))
        except (ValueError, TypeError):
            return "- &deg;C"

    context = {
        'name': result.name or '-',
        'chromosome': result.chromosome,
        'start': result.start,
        'end': result.end,
        'hrfw': result.hr_dna[0],
        'hrrv': result.hr_dna[1],
        'deleted_dna': result.hr_dna[2],
        'primer_left': primer.get('PRIMER_LEFT_0_SEQUENCE', '-'),
        'left_tm': get_tm_str(primer.get('PRIMER_LEFT_0_TM', 0)),
        'primer_right': primer.get('PRIMER_RIGHT_0_SEQUENCE', '-'),
        'right_tm': get_tm_str(primer.get('PRIMER_RIGHT_0_TM', 0)),
        'deleted_dna_size': str(
            primer.get('PRIMER_PAIR_0_PRODUCT_SIZE', '-')
        ) + " (bp)",
        'negative_result_size': str(
            primer.get('negative_result', '-')
        ) + " (bp)",
    }
    context['json_table'] = json.dumps(result.guide_table)
    context['spedit_json'] = json.dumps(
        build_spedit_candidate_data(result.guide_table)
    )

    return template_text % context
