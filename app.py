import csv
import operator
from flask import Flask
from flask_restplus import Resource, Api, apidoc
from flask import jsonify
from werkzeug.contrib.fixers import ProxyFix
from rdkit import Chem
from rdkit_functions import get_rdkit_descriptors
from similarity_functions import similar_FDA_molecules
#from rdkit import DataStructs
from rdkit.Chem import Descriptors
#from rdkit.Chem.Fingerprints import FingerprintMols
from sanifix import fix_mol

# Authentication
from functools import wraps
from flask import request, Response

def check_auth(username, password):
    return username == 'demo' and password == 'demo1234'

def authenticate():
    message = {'message': "Could not verify your access level for that URL.\nPlease login with proper credentials."}
    resp = jsonify(message)

    resp.status_code = 401
    resp.headers['WWW-Authenticate'] = 'Basic realm="Example"'

    return resp

def requires_auth(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        auth = request.authorization
        if not auth: 
            return authenticate()

        elif not check_auth(auth.username, auth.password):
            return authenticate()
        return f(*args, **kwargs)

    return decorated


app = Flask(__name__)
app.config.SWAGGER_UI_DOC_EXPANSION = 'list'
api = Api(app, version='0.1', title='Explorable Labs API', description='API microservices for science.\nThe molecular descriptors and some of the similarity-related operations are being generated using RDKit version 2018.09.1\nThe json returned from these endpoints is meant to be consumed by code - if you want it to be human-readable you should install an extension on your browser for viewing JSON files, such as JSON View.')
api.namespaces = []
app.wsgi_app = ProxyFix(app.wsgi_app) # Let the swagger docs show up when we're on https

ns = api.namespace('descriptors', description="Operations related to molecular descriptors. Example <a href=\"https://api.explorablelabs.com/descriptors/smiles/CN1C=NC2=C1C(=O)N(C(=O)N2C)C\">caffeine</a>.")
ns_sim = api.namespace('similarity', description="Operations related to molecular similarity. (Requires authentication. Contact support@explorablelabs.com for a demo login.) Example <a href=\"https://api.explorablelabs.com/similarity/FDA/smiles/CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5\">Gleevec</a>.")

@ns.route('/smiles/<string:smiles>')
@api.response(404, 'Not found')
class DescribeSmiles(Resource):
    def get(self, smiles):
        """Returns molecular descriptors for a valid SMILES string"""
        m = Chem.MolFromSmiles(smiles)
        m = fix_mol(m)
        return jsonify(get_rdkit_descriptors(m))

@ns.route('/', methods = ['POST'])
def get_smiles_descriptors():
    if request.headers['Content-Type'] == 'text/plain':
        smiles = request.data
        m = None
        try:
            m = Chem.MolFromSmiles(smiles)
        except:
            m = None
        if not m:
            try:
                m = Chem.MolFromSmarts(smarts)
            except:
                m = None
        if not m:
            return jsonify("Not a valid SMILES string.")
        else:
            m = fix_mol(m)
            return jsonify(get_rdkit_descriptors(m))
    else:
        return jsonify("Not a valid SMILES string.")

@ns.route('/smarts/<string:smarts>')
@api.response(404, 'Not found')
class DescribeSmarts(Resource):
    def get(self, smarts):
        """Returns molecular descriptors for a valid SMARTS string"""
        m = Chem.MolFromSmarts(smarts)
        m = fix_mol(m)
        return jsonify(get_rdkit_descriptors(m))

@ns_sim.route('/FDA/smiles/<string:smiles>')
@api.response(404, 'Not found')
class SimilarFDA(Resource):
    @requires_auth
    def get(self, smiles):
        """Returns FDA compounds most similar to a valid SMILES string"""
        m = Chem.MolFromSmiles(smiles)
        m = fix_mol(m)
        return jsonify(similar_FDA_molecules(m))

@ns_sim.route('/FDA/smarts/<string:smarts>')
@api.response(404, 'Not found')
class SimilarFDA(Resource):
    def get(self, smarts):
        """Returns FDA compounds most similar to a valid SMARTS string"""
        m = Chem.MolFromSmarts(smarts)
        m = fix_mol(m)
        return jsonify(similar_FDA_molecules(m))

if __name__ == '__main__':
    app.run(debug=False, host="0.0.0.0")

