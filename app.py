from flask import Flask, render_template, request, send_file, flash, redirect, url_for, send_from_directory, render_template_string
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
from werkzeug.utils import secure_filename
from PIL import Image, ImageDraw, ImageFont

app = Flask(__name__)
app.secret_key = 'your_secret_key_here'

# Ensure the upload and results folders exist
UPLOAD_FOLDER = 'uploads'
RESULTS_FOLDER = 'results'
STATIC_FOLDER = 'static'
for folder in [UPLOAD_FOLDER, RESULTS_FOLDER, STATIC_FOLDER]:
    if not os.path.exists(folder):
        os.makedirs(folder)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['RESULTS_FOLDER'] = RESULTS_FOLDER
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16MB max-limit

ALLOWED_EXTENSIONS = {'ab1', 'fasta', 'fa'}

IUPAC_CODES = {
    "R": ["A", "G"], "Y": ["C", "T"], "S": ["G", "C"], "W": ["A", "T"],
    "K": ["G", "T"], "M": ["A", "C"], "B": ["C", "G", "T"], "D": ["A", "G", "T"],
    "H": ["A", "C", "T"], "V": ["A", "C", "G"], "N": ["A", "C", "G", "T"]
}

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def create_text_favicon(text, size=(32, 32), bg_color=(0, 0, 0), text_color=(255, 255, 255)):
    image = Image.new('RGB', size, bg_color)
    draw = ImageDraw.Draw(image)
    try:
        font = ImageFont.truetype("arial.ttf", 12)
    except IOError:
        font = ImageFont.load_default()
    text_width, text_height = draw.textsize(text, font=font)
    position = ((size[0] - text_width) / 2, (size[1] - text_height) / 2)
    draw.text(position, text, font=font, fill=text_color)
    favicon_path = os.path.join(STATIC_FOLDER, 'favicon.ico')
    image.save(favicon_path)
    return favicon_path

# Create favicon if it doesn't exist
if not os.path.exists(os.path.join(STATIC_FOLDER, 'favicon.ico')):
    create_text_favicon("MSA")

def format_alignment(ref_seq, ab1_seq, alignment):
    aligned_ref = alignment[0]
    aligned_ab1 = alignment[1]
    
    result = ""
    line_width = 60
    
    for i in range(0, len(aligned_ref), line_width):
        ref_line = aligned_ref[i:i+line_width]
        ab1_line = aligned_ab1[i:i+line_width]
        
        # Create conservation line
        conservation = ''.join(['*' if r == a and r != '-' else ' ' for r, a in zip(ref_line, ab1_line)])
        
        result += f"Ref sequence  - {ref_line}\n"
        result += f"Ab1 Sequence  - {ab1_line}\n"
        result += f"Conservation  - {conservation}\n\n"
    
    return result.strip()

def resolve_ambiguity(base1, base2):
    """Check if two bases match, including handling IUPAC ambiguous codes."""
    if base1 == base2:
        return True
    if base1 in IUPAC_CODES and base2 in IUPAC_CODES[base1]:
        return True
    if base2 in IUPAC_CODES and base1 in IUPAC_CODES[base2]:
        return True
    return False

@app.route('/', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        # Check if both files are present in the request
        if 'ab1_file' not in request.files or 'reference_file' not in request.files:
            flash('Both AB1 and reference files are required')
            return redirect(request.url)
        
        ab1_file = request.files['ab1_file']
        reference_file = request.files['reference_file']
        
        # Check if filenames are empty
        if ab1_file.filename == '' or reference_file.filename == '':
            flash('No selected file')
            return redirect(request.url)
        
        if ab1_file and reference_file and allowed_file(ab1_file.filename) and allowed_file(reference_file.filename):
            # Save uploaded files
            ab1_filename = secure_filename(ab1_file.filename)
            ref_filename = secure_filename(reference_file.filename)
            ab1_path = os.path.join(app.config['UPLOAD_FOLDER'], ab1_filename)
            ref_path = os.path.join(app.config['UPLOAD_FOLDER'], ref_filename)
            ab1_file.save(ab1_path)
            reference_file.save(ref_path)
            
            try:
                # Process the AB1 file
                ab1_record = SeqIO.read(ab1_path, "abi")
                ab1_seq = str(ab1_record.seq)
                
                # Process the reference file
                ref_record = SeqIO.read(ref_path, "fasta")
                ref_seq = str(ref_record.seq)
                
                # Perform pairwise alignment
                alignments = pairwise2.align.globalcs(
                    ref_seq, ab1_seq,
                    lambda x, y: 2 if resolve_ambiguity(x, y) else -1,  # Match score if bases align
                    -0.5,  # Gap opening penalty
                    -0.1   # Gap extension penalty
                )
                best_alignment = alignments[0]
                
                # Format the alignment
                alignment_result = format_alignment(ref_seq, ab1_seq, best_alignment)
                
                # Save the alignment result
                result_filename = 'alignment_result.txt'
                result_path = os.path.join(app.config['RESULTS_FOLDER'], result_filename)
                with open(result_path, 'w') as f:
                    f.write(alignment_result)
                
                # Redirect to the result page with the alignment result
                return redirect(url_for('result', alignment_result=alignment_result))
            
            except Exception as e:
                flash(f'An error occurred: {str(e)}')
                return redirect(request.url)
        else:
            flash('Invalid file type. Allowed file types are: .ab1, .fasta, .fa')
            return redirect(request.url)
    
    return render_template('upload.html')

@app.route('/favicon.ico')
def favicon():
    return send_from_directory(os.path.join(app.root_path, 'static'),
                               'favicon.ico', mimetype='image/vnd.microsoft.icon')

# New route for displaying results
@app.route('/result')
def result():
    alignment_result = request.args.get('alignment_result', None)  # Get the alignment result from the query parameter
    return render_template('result.html', alignment_result=alignment_result)  # Pass the result to the template

if __name__ == '__main__':
    app.run(debug=True)
