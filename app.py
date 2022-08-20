from flask import Flask, render_template, request, redirect, jsonify
from flask_wtf import FlaskForm
from wtforms import StringField, FileField, SubmitField, RadioField, TextAreaField
from wtforms.validators import DataRequired
from Blastlogica import *
import redis

app = Flask(__name__)
cache = redis.Redis(host='redis', port=6379)
app.config['SECRET_KEY'] = 'mykey'
app.config["UPLOAD_FOLDER"] = "tmp/"

#genereert user input form objecten voor de HTML pagina.
class InputForm(FlaskForm):
    sequence = TextAreaField()
    file = FileField('File')
    seqType = RadioField('Sequentietype:',
                            choices=[('nuc', 'DNA/RNA'), ('prot', 'Proteine')],validators=[DataRequired()])
    ziekte = StringField("Ziekte:")
    submit = SubmitField('Submit')

#REST API methode; accepteert json/dictionary van 1 of meerdere sequenties met sequentietype ("nuc", "prot")
#en evt een ziekte als value, retourneert een jsonfile van de annotatie.
#redirects naar de webpagina /Blasttool wanneer er geen jsonfile aangeleverd wordt.
@app.route('/', methods = ['GET', 'POST'])
def main():
    if request.method == 'GET':
        return redirect('/BlastTool')
    if request.method == 'POST':
        i = 1
        blast_records = {}
        response = request.get_json()
        for sequence in response:
            blast_records[i] = seqBlastRecords(blast_seq(sequence, response[sequence]["seqType"]))
            i += 1
        return jsonify(blast_records)

#methode voor renderen van de webpaginas en het accepteren en verwerken en retourneren van user input.
@app.route('/BlastTool', methods = ['GET', 'POST'])
def index():
    error = None
    form = InputForm(request.form)

    if form.validate_on_submit():
        input_seq = form.sequence.data
        seqType = form.seqType.data
        ziekte = form.ziekte.data
        file = request.files['file']

        if not file and input_seq == "":
            error = "Voer een sequentie of fastfile in"
            return render_template("index.html", form=form, error=error)
        elif input_seq != "":
            i = 1
            blast_records = {}
            sequences = input_seq.split("\n")
            if "" in sequences:
                sequences.remove("")
            for seq in sequences:
                blast_records[i] = seqBlastRecords(blast_seq(seq, seqType))
                i += 1
            return render_template('content.html', blast_records=blast_records)
        else:
            file.save(app.config['UPLOAD_FOLDER']+ "fasta_upload")
            blast_records = fastaBlastRecords(blast_fasta(seqType))
            return render_template('content.html', blast_records=blast_records)

    return render_template('index.html', form=form)

#methode die het aangeleverde .fasta bestand inleest, met qBlast runt en parsed.
#retourneert een NCBIXML parsed bestand van de alignments.
def blast_fasta(seq_type):
    fasta_string = open(app.config['UPLOAD_FOLDER'] + "fasta_upload","r").read()
    if seq_type == "nuc":
        result_handle = NCBIXML.parse(NCBIWWW.qblast("blastn", "nt", fasta_string))
    else:
        result_handle = NCBIXML.parse(NCBIWWW.qblast("blastp", "nr", fasta_string))
    return result_handle

if __name__ == '__main__':
    app.run(debug = True)
