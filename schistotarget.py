from flask import Flask,flash,render_template,request,redirect, url_for
from werkzeug.contrib.fixers import ProxyFix
import os, sys, subprocess
import time
import sqlite3
import hashlib
from werkzeug import secure_filename
sys.path.append(os.getcwd()+'/scripts')
from features import features
from Bio import SeqIO
from StringIO import StringIO
from celery import Celery
from celery.result import AsyncResult
from sklearn import svm, preprocessing
import numpy as np


#----------------------CELERY CONFIGURATION-----------------------

app = Flask(__name__)
app.config['CELERY_ACCEPT_CONTENT'] = ['json','pickle']
app.config['CELERY_TASK_SERIALIZER'] = 'json'
app.config['CELERY_BROKER_URL'] = 'redis://localhost:6379/0'
app.config['CELERY_RESULT_BACKEND'] = 'redis://localhost:6379/0'

celery = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])
celery.conf.update(app.config)


#----------------------IMMUNO-REACTIVE PROTEINS PEDICTION-----------------------
@celery.task()
def run_immuno(filename, immuno_email):
    conn = sqlite3.connect('database.db')
    c = conn.cursor()
    
    parameter=""
    seqID_list=[]
    result_file=open(filename+"_result.txt", 'w')
    result_file.write("Sequence_ID\tPrediction\tProbability\n")

    records=SeqIO.parse(filename, "fasta")
    for record in records:
        hash_sequence=hashlib.md5(str(record.seq)).hexdigest()
        c.execute("SELECT * FROM immuno WHERE sequence='"+hash_sequence+"'")
        data=c.fetchone()
        if data is None:
            parameter=parameter+features(record.id, str(record.seq))+"\n"
            seqID_list.append(record.id)
                 
        else:
            c.execute("UPDATE immuno SET access=access+1, time=CURRENT_TIMESTAMP WHERE sequence='"+hash_sequence+"'")
            conn.commit()
            c.execute("SELECT prediction, score FROM immuno WHERE sequence='"+hash_sequence+"'")
            data1=c.fetchone()
            result_file.write(str(record.id)+"\t"+data1[0]+"\t"+str(data1[1])+"\n")


    records.close()
    
    
    #---------------------WORKING WITH SCIKIT-LEARN FOR IMMUNO-REACTIVE PROTEINS------------
    if parameter!="":
        parameter=StringIO(parameter)
        train_para="immuno.csv"
        train_label="immuno_labels.csv"
        
        train_data=np.genfromtxt(train_para, delimiter=",")
        

        train_label=np.genfromtxt(train_label, delimiter=",")
        test_data=np.genfromtxt(parameter, delimiter=",")
        
        min_max_scaler = preprocessing.MinMaxScaler(feature_range=(0, 1))

        train_data_scaled = min_max_scaler.fit_transform(train_data)

        test_data_scaled = min_max_scaler.fit_transform(test_data)

        clf = svm.SVC(kernel='rbf', C=10.0, gamma=2.0, probability=True)

        clf.fit(train_data_scaled, train_label)

        test_predicted = clf.predict(test_data_scaled)

        scores = clf.predict_proba(test_data_scaled)
        

        #-----------------
        fasta_rec=SeqIO.index(filename, "fasta")
        i=0
        for pred in test_predicted:
            score=scores[:,1][i]
            if pred==1.0:
                pred='Immunoreactivity'
            if pred==0.0:
                pred='No Immunoreactivity'

            if score>0.5:
                result_file.write(seqID_list[i]+"\t"+pred+"\t"+str(score)+"\n")
                c.execute("INSERT INTO immuno VALUES ('"+hashlib.md5(str(fasta_rec[seqID_list[i]].seq)).hexdigest()+"', '"+pred+"', '"+str(score)+"', 0, CURRENT_TIMESTAMP)")
          
            else:
                result_file.write(seqID_list[i]+"\t"+pred+"\t"+str(1-score)+"\n")
                c.execute("INSERT INTO immuno VALUES ('"+hashlib.md5(str(fasta_rec[seqID_list[i]].seq)).hexdigest()+"', '"+pred+"', '"+str(1-score)+"', 0, CURRENT_TIMESTAMP)")
            
            i=i+1
        conn.commit()
        fasta_rec.close()
    result_file.close()

    #--------------------EMAIL SENDING--------------------------------------

    if immuno_email!="":
        command = "echo 'Your SchistoTarget Prediction Result is attached for job ID: "+filename+"\n\n\nKind regards,\n\nLutz Krause & Shihab Hasan\nComputational Medical Genomics Group, The University of Queensland Diamantina Institute' | mutt -a "+filename+"'_result.txt' -s 'SchistoTarget Prediction Result' -- "+immuno_email
        subprocess.call(command, shell=(sys.platform!="Linux"))

    f = open(filename+"_result.txt",'r')
    lines = f.readlines()[1:]
    f.close()
    os.remove(filename)
    os.remove(filename+"_result.txt")
    return lines

#----------------------IgE ProteinS PEDICTION-----------------------
@celery.task()
def run_IgE(filename, IgE_email):
    conn = sqlite3.connect('database.db')
    c = conn.cursor()
    
    parameter=""
    seqID_list=[]
    result_file=open(filename+"_result.txt", 'w')
    result_file.write("Sequence_ID\tPrediction\tProbability\n")

    records=SeqIO.parse(filename, "fasta")
    for record in records:
        hash_sequence=hashlib.md5(str(record.seq)).hexdigest()
        c.execute("SELECT * FROM IgE WHERE sequence='"+hash_sequence+"'")
        data=c.fetchone()
        if data is None:
            parameter=parameter+features(record.id, str(record.seq))+"\n"
            seqID_list.append(record.id)
                 
        else:
            c.execute("UPDATE IgE SET access=access+1, time=CURRENT_TIMESTAMP WHERE sequence='"+hash_sequence+"'")
            conn.commit()
            c.execute("SELECT prediction, score FROM IgE WHERE sequence='"+hash_sequence+"'")
            data1=c.fetchone()
            result_file.write(str(record.id)+"\t"+data1[0]+"\t"+str(data1[1])+"\n")


    records.close()
    
    
    #---------------------WORKING WITH SCIKIT-LEARN FOR IgE------------
    if parameter!="":
        parameter=StringIO(parameter)
        train_para="IgE.csv"
        train_label="IgE_labels.csv"
        
        train_data=np.genfromtxt(train_para, delimiter=",")
        

        train_label=np.genfromtxt(train_label, delimiter=",")
        test_data=np.genfromtxt(parameter, delimiter=",")
        
        min_max_scaler = preprocessing.MinMaxScaler(feature_range=(0, 1))

        train_data_scaled = min_max_scaler.fit_transform(train_data)

        test_data_scaled = min_max_scaler.fit_transform(test_data)

        clf = svm.SVC(kernel='rbf', C=8.0, gamma=2.0, probability=True)

        clf.fit(train_data_scaled, train_label)

        test_predicted = clf.predict(test_data_scaled)

        scores = clf.predict_proba(test_data_scaled)
        

        #-----------------
        fasta_rec=SeqIO.index(filename, "fasta")
        i=0
        for pred in test_predicted:
            score=scores[:,1][i]
            if pred==1.0:
                pred='IgE reactivity'
            if pred==0.0:
                pred='No IgE reactivity'

            if score>0.5:
                result_file.write(seqID_list[i]+"\t"+pred+"\t"+str(score)+"\n")
                c.execute("INSERT INTO IgE VALUES ('"+hashlib.md5(str(fasta_rec[seqID_list[i]].seq)).hexdigest()+"', '"+pred+"', '"+str(score)+"', 0, CURRENT_TIMESTAMP)")
          
            else:
                result_file.write(seqID_list[i]+"\t"+pred+"\t"+str(1-score)+"\n")
                c.execute("INSERT INTO IgE VALUES ('"+hashlib.md5(str(fasta_rec[seqID_list[i]].seq)).hexdigest()+"', '"+pred+"', '"+str(1-score)+"', 0, CURRENT_TIMESTAMP)")
            
            i=i+1
        conn.commit()
        fasta_rec.close()
    result_file.close()

    if IgE_email!="":
        command = "echo 'Your SchistoTarget Prediction Result is attached for job ID: '"+filename+"'\n\n\nKind regards,\n\nLutz Krause & Shihab Hasan\nComputational Medical Genomics Group, The University of Queensland Diamantina Institute'"+" | mutt -a "+filename+"'_result.txt' -s 'SchistoTarget Prediction Result' -- "+IgE_email
        subprocess.call(command, shell=(sys.platform!="Linux"))

    f = open(filename+"_result.txt",'r')
    lines = f.readlines()[1:]
    f.close()
    os.remove(filename)
    os.remove(filename+"_result.txt")
    return lines


	
#----------------------IgG1 ProteinS PEDICTION-----------------------
@celery.task()
def run_IgG1(filename, IgG1_email):
    conn = sqlite3.connect('database.db')
    c = conn.cursor()
    
    parameter=""
    seqID_list=[]
    result_file=open(filename+"_result.txt", 'w')
    result_file.write("Sequence_ID\tPrediction\tProbability\n")

    records=SeqIO.parse(filename, "fasta")
    for record in records:
        hash_sequence=hashlib.md5(str(record.seq)).hexdigest()
        c.execute("SELECT * FROM IgG1 WHERE sequence='"+hash_sequence+"'")
        data=c.fetchone()
        if data is None:
            parameter=parameter+features(record.id, str(record.seq))+"\n"
            seqID_list.append(record.id)
                 
        else:
            c.execute("UPDATE IgG1 SET access=access+1, time=CURRENT_TIMESTAMP WHERE sequence='"+hash_sequence+"'")
            conn.commit()
            c.execute("SELECT prediction, score FROM IgG1 WHERE sequence='"+hash_sequence+"'")
            data1=c.fetchone()
            result_file.write(str(record.id)+"\t"+data1[0]+"\t"+str(data1[1])+"\n")


    records.close()
    
    
    #---------------------WORKING WITH SCIKIT-LEARN FOR IgG1------------
    if parameter!="":
        parameter=StringIO(parameter)
        train_para="IgG1.csv"
        train_label="IgG1_labels.csv"
        
        train_data=np.genfromtxt(train_para, delimiter=",")
        

        train_label=np.genfromtxt(train_label, delimiter=",")
        test_data=np.genfromtxt(parameter, delimiter=",")
        
        min_max_scaler = preprocessing.MinMaxScaler(feature_range=(0, 1))

        train_data_scaled = min_max_scaler.fit_transform(train_data)

        test_data_scaled = min_max_scaler.fit_transform(test_data)

        clf = svm.SVC(kernel='rbf', C=8.0, gamma=2.0, probability=True)

        clf.fit(train_data_scaled, train_label)

        test_predicted = clf.predict(test_data_scaled)

        scores = clf.predict_proba(test_data_scaled)
        

        #-----------------
        fasta_rec=SeqIO.index(filename, "fasta")
        i=0
        for pred in test_predicted:
            score=scores[:,1][i]
            if pred==1.0:
                pred='IgG1 reactivity'
            if pred==0.0:
                pred='No IgG1 reactivity'

            if score>0.5:
                result_file.write(seqID_list[i]+"\t"+pred+"\t"+str(score)+"\n")
                c.execute("INSERT INTO IgG1 VALUES ('"+hashlib.md5(str(fasta_rec[seqID_list[i]].seq)).hexdigest()+"', '"+pred+"', '"+str(score)+"', 0, CURRENT_TIMESTAMP)")
          
            else:
                result_file.write(seqID_list[i]+"\t"+pred+"\t"+str(1-score)+"\n")
                c.execute("INSERT INTO IgG1 VALUES ('"+hashlib.md5(str(fasta_rec[seqID_list[i]].seq)).hexdigest()+"', '"+pred+"', '"+str(1-score)+"', 0, CURRENT_TIMESTAMP)")
            
            i=i+1
        conn.commit()
        fasta_rec.close()
    result_file.close()

    if IgG1_email!="":
        command = "echo 'Your SchistoTarget Prediction Result is attached for job ID: '"+filename+"'\n\n\nKind regards,\n\nLutz Krause & Shihab Hasan\nComputational Medical Genomics Group, The University of Queensland Diamantina Institute'"+" | mutt -a "+filename+"'_result.txt' -s 'SchistoTarget Prediction Result' -- "+IgG1_email
        subprocess.call(command, shell=(sys.platform!="Linux"))

    f = open(filename+"_result.txt",'r')
    lines = f.readlines()[1:]
    f.close()
    os.remove(filename)
    os.remove(filename+"_result.txt")
    return lines





#----------------------IgG3 ProteinS PEDICTION-----------------------
@celery.task()
def run_IgG3(filename, IgG3_email):
    conn = sqlite3.connect('database.db')
    c = conn.cursor()
    
    parameter=""
    seqID_list=[]
    result_file=open(filename+"_result.txt", 'w')
    result_file.write("Sequence_ID\tPrediction\tProbability\n")

    records=SeqIO.parse(filename, "fasta")
    for record in records:
        hash_sequence=hashlib.md5(str(record.seq)).hexdigest()
        c.execute("SELECT * FROM IgG3 WHERE sequence='"+hash_sequence+"'")
        data=c.fetchone()
        if data is None:
            parameter=parameter+features(record.id, str(record.seq))+"\n"
            seqID_list.append(record.id)
                 
        else:
            c.execute("UPDATE IgG3 SET access=access+1, time=CURRENT_TIMESTAMP WHERE sequence='"+hash_sequence+"'")
            conn.commit()
            c.execute("SELECT prediction, score FROM IgG3 WHERE sequence='"+hash_sequence+"'")
            data1=c.fetchone()
            result_file.write(str(record.id)+"\t"+data1[0]+"\t"+str(data1[1])+"\n")


    records.close()
    
    
    #---------------------WORKING WITH SCIKIT-LEARN FOR IgG3------------
    if parameter!="":
        parameter=StringIO(parameter)
        train_para="IgG3.csv"
        train_label="IgG3_labels.csv"
        
        train_data=np.genfromtxt(train_para, delimiter=",")
        

        train_label=np.genfromtxt(train_label, delimiter=",")
        test_data=np.genfromtxt(parameter, delimiter=",")
        
        min_max_scaler = preprocessing.MinMaxScaler(feature_range=(0, 1))

        train_data_scaled = min_max_scaler.fit_transform(train_data)

        test_data_scaled = min_max_scaler.fit_transform(test_data)

        clf = svm.SVC(kernel='rbf', C=10.0, gamma=10.0, probability=True)

        clf.fit(train_data_scaled, train_label)

        test_predicted = clf.predict(test_data_scaled)

        scores = clf.predict_proba(test_data_scaled)
        

        #-----------------
        fasta_rec=SeqIO.index(filename, "fasta")
        i=0
        for pred in test_predicted:
            score=scores[:,1][i]
            if pred==1.0:
                pred='IgG3 reactivity'
            if pred==0.0:
                pred='No IgG3 reactivity'

            if score>0.5:
                result_file.write(seqID_list[i]+"\t"+pred+"\t"+str(score)+"\n")
                c.execute("INSERT INTO IgG3 VALUES ('"+hashlib.md5(str(fasta_rec[seqID_list[i]].seq)).hexdigest()+"', '"+pred+"', '"+str(score)+"', 0, CURRENT_TIMESTAMP)")
          
            else:
                result_file.write(seqID_list[i]+"\t"+pred+"\t"+str(1-score)+"\n")
                c.execute("INSERT INTO IgG3 VALUES ('"+hashlib.md5(str(fasta_rec[seqID_list[i]].seq)).hexdigest()+"', '"+pred+"', '"+str(1-score)+"', 0, CURRENT_TIMESTAMP)")
            
            i=i+1
        conn.commit()
        fasta_rec.close()
    result_file.close()

    if IgG3_email!="":
        command = "echo 'Your SchistoTarget Prediction Result is attached for job ID: '"+filename+"'\n\n\nKind regards,\n\nLutz Krause & Shihab Hasan\nComputational Medical Genomics Group, The University of Queensland Diamantina Institute'"+" | mutt -a "+filename+"'_result.txt' -s 'SchistoTarget Prediction Result' -- "+IgG3_email
        subprocess.call(command, shell=(sys.platform!="Linux"))

    f = open(filename+"_result.txt",'r')
    lines = f.readlines()[1:]
    f.close()
    os.remove(filename)
    os.remove(filename+"_result.txt")
    return lines




#----------------------IgG4 ProteinS PEDICTION-----------------------
@celery.task()
def run_IgG4(filename, IgG4_email):
    conn = sqlite3.connect('database.db')
    c = conn.cursor()
    
    parameter=""
    seqID_list=[]
    result_file=open(filename+"_result.txt", 'w')
    result_file.write("Sequence_ID\tPrediction\tProbability\n")

    records=SeqIO.parse(filename, "fasta")
    for record in records:
        hash_sequence=hashlib.md5(str(record.seq)).hexdigest()
        c.execute("SELECT * FROM IgG4 WHERE sequence='"+hash_sequence+"'")
        data=c.fetchone()
        if data is None:
            parameter=parameter+features(record.id, str(record.seq))+"\n"
            seqID_list.append(record.id)
                 
        else:
            c.execute("UPDATE IgG4 SET access=access+1, time=CURRENT_TIMESTAMP WHERE sequence='"+hash_sequence+"'")
            conn.commit()
            c.execute("SELECT prediction, score FROM IgG4 WHERE sequence='"+hash_sequence+"'")
            data1=c.fetchone()
            result_file.write(str(record.id)+"\t"+data1[0]+"\t"+str(data1[1])+"\n")


    records.close()
    
    
    #---------------------WORKING WITH SCIKIT-LEARN FOR IgG4------------
    if parameter!="":
        parameter=StringIO(parameter)
        train_para="IgG4.csv"
        train_label="IgG4_labels.csv"
        
        train_data=np.genfromtxt(train_para, delimiter=",")
        

        train_label=np.genfromtxt(train_label, delimiter=",")
        test_data=np.genfromtxt(parameter, delimiter=",")
        
        min_max_scaler = preprocessing.MinMaxScaler(feature_range=(0, 1))

        train_data_scaled = min_max_scaler.fit_transform(train_data)

        test_data_scaled = min_max_scaler.fit_transform(test_data)

        clf = svm.SVC(kernel='rbf', C=1.0, gamma=1.0, probability=True)

        clf.fit(train_data_scaled, train_label)

        test_predicted = clf.predict(test_data_scaled)

        scores = clf.predict_proba(test_data_scaled)
        

        #-----------------
        fasta_rec=SeqIO.index(filename, "fasta")
        i=0
        for pred in test_predicted:
            score=scores[:,1][i]
            if pred==1.0:
                pred='IgG4 reactivity'
            if pred==0.0:
                pred='No IgG4 reactivity'

            if score>0.5:
                result_file.write(seqID_list[i]+"\t"+pred+"\t"+str(score)+"\n")
                c.execute("INSERT INTO IgG4 VALUES ('"+hashlib.md5(str(fasta_rec[seqID_list[i]].seq)).hexdigest()+"', '"+pred+"', '"+str(score)+"', 0, CURRENT_TIMESTAMP)")
          
            else:
                result_file.write(seqID_list[i]+"\t"+pred+"\t"+str(1-score)+"\n")
                c.execute("INSERT INTO IgG4 VALUES ('"+hashlib.md5(str(fasta_rec[seqID_list[i]].seq)).hexdigest()+"', '"+pred+"', '"+str(1-score)+"', 0, CURRENT_TIMESTAMP)")
            
            i=i+1
        conn.commit()
        fasta_rec.close()
    result_file.close()

    if IgG4_email!="":
        command = "echo 'Your SchistoTarget Prediction Result is attached for job ID: '"+filename+"'\n\n\nKind regards,\n\nLutz Krause & Shihab Hasan\nComputational Medical Genomics Group, The University of Queensland Diamantina Institute'"+" | mutt -a "+filename+"'_result.txt' -s 'SchistoTarget Prediction Result' -- "+IgG4_email
        subprocess.call(command, shell=(sys.platform!="Linux"))

    f = open(filename+"_result.txt",'r')
    lines = f.readlines()[1:]
    f.close()
    os.remove(filename)
    os.remove(filename+"_result.txt")
    return lines




#----------------------WORKING WITH HTML-----------------------------------
app = Flask(__name__)
UPLOAD_FOLDER = os.getcwd()
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.secret_key = 'some_secret'

@app.route('/',methods=['GET','POST'])
def index():
   return render_template('index.html')

@app.route('/home',methods=['GET','POST'])
def home():
   return render_template('index.html')

@app.route('/help',methods=['GET','POST'])
def manual():
   return render_template('help.html')

@app.route('/contact',methods=['GET','POST'])
def contact():
    return render_template('contact.html')

@app.route('/thanks',methods=['POST'])
def thanks():
    name=request.form['name']
    email=request.form['email']
    message=request.form['message']
    command = "echo 'Name: '"+name+"'\nEmail: '"+email+"'\nMessage: '"+message+" | mutt -s 'SchistoTarget Prediction Query' -- pharm.shihab@gmail.com"
    subprocess.call(command, shell=(sys.platform!="Linux"))
    return render_template('thanks.html', name=name)

@app.route('/immunoreactivity_app',methods=['GET','POST'])
def immunoreactivity_app():
    return render_template('immunoreactivity.html')
	
@app.route('/IgE_app',methods=['GET','POST'])
def IgE_app():
    return render_template('IgE.html')
	
@app.route('/IgG1_app',methods=['GET','POST'])	
def IgG1_app():
    return render_template('IgG1.html')

@app.route('/IgG3_app',methods=['GET','POST'])
def IgG3_app():
    return render_template('IgG3.html')
	
@app.route('/IgG4_app',methods=['GET','POST'])
def IgG4_app():
    return render_template('IgG4.html')

#----------------------immuno-------------------------

@app.route('/immuno_predict',methods=['POST'])
def immuno_predict():
    global filename
    global immuno_email
    immuno_email=request.form['immuno_email'].replace(" ","")
    if request.form['immuno_sequences'].replace(" ","")!="":
        filename=hashlib.md5(time.asctime()).hexdigest()
        file_in=open(filename, 'w')
        file_in.write(request.form['immuno_sequences'].replace(" >",">"))
        file_in.close()
    else:
        file = request.files['immuno_file']
        filename = secure_filename(file.filename)+"_"+hashlib.md5(time.asctime()).hexdigest()
        file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))

    task=run_immuno.delay(filename, immuno_email)
    id=task.id
    return redirect(url_for('immuno_progress', id=id))



@app.route('/immuno_progress/<id>')
def immuno_progress(id):
    return render_template('immuno_progress.html', id=id)
      
@app.route('/immuno_status/<id>')
def immuno_status(id):
    """ Wait until the job is finished and report success."""
    global lines
    result = run_immuno.AsyncResult(id, app=celery)
    lines=result.get()
    while not result.ready():
        if result.failed():
            flash(u"Error running SchistoTarget. please check input file", 'error')
            return redirect(url_for('index'))
        time.sleep(5)
    return 'success'

@app.route('/immuno_results/<id>')
def immuno_results(id):
    return render_template('immuno_result.html', immuno_result=lines)

#----------------------IgE-------------------------

@app.route('/IgE_predict',methods=['POST'])
def IgE_predict():   
    global filename
    global IgE_email
    IgE_email=request.form['IgE_email'].replace(" ","")
    if request.form['IgE_sequences'].replace(" ","")!="":
        filename=hashlib.md5(time.asctime()).hexdigest()
        file_in=open(filename, 'w')
        file_in.write(request.form['IgE_sequences'].replace(" >",">"))
        file_in.close()
    else:
        file = request.files['IgE_file']
        filename = secure_filename(file.filename)+"_"+hashlib.md5(time.asctime()).hexdigest()
        file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))

    task=run_IgE.delay(filename, IgE_email)
    id=task.id
    return redirect(url_for('IgE_progress', id=id))

@app.route('/IgE_progress/<id>')
def IgE_progress(id):
    return render_template('IgE_progress.html', id=id)

@app.route('/IgE_status/<id>')
def IgE_status(id):
    """ Wait until the job is finished and report success."""
    global lines
    result = run_IgE.AsyncResult(id, app=celery)
    lines=result.get()
    while not result.ready():
        if result.failed():
            flash(u"Error running SchistoTarget. please check input file", 'error')
            return redirect(url_for('index'))
        time.sleep(5)
    return 'success'

@app.route('/IgE_results/<id>')
def IgE_results(id):
    return render_template('IgE_result.html', IgE_result=lines)
	
	
#----------------------IgG1-------------------------
@app.route('/IgG1_predict',methods=['POST'])
def IgG1_predict():
    global filename
    global IgG1_email
    IgG1_email=request.form['IgG1_email'].replace(" ","")
    if request.form['IgG1_sequences'].replace(" ","")!="":
        filename=hashlib.md5(time.asctime()).hexdigest()
        file_in=open(filename, 'w')
        file_in.write(request.form['IgG1_sequences'].replace(" >",">"))
        file_in.close()
    else:
        file = request.files['IgG1_file']
        filename = secure_filename(file.filename)+"_"+hashlib.md5(time.asctime()).hexdigest()
        file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))

    task=run_IgG1.delay(filename, IgG1_email)
    id=task.id
    return redirect(url_for('IgG1_progress', id=id))



@app.route('/IgG1_progress/<id>')
def IgG1_progress(id):
    return render_template('IgG1_progress.html', id=id)
      
@app.route('/IgG1_status/<id>')
def IgG1_status(id):
    """ Wait until the job is finished and report success."""
    global lines
    result = run_IgG1.AsyncResult(id, app=celery)
    lines=result.get()
    while not result.ready():
        if result.failed():
            flash(u"Error running SchistoProt. please check input file", 'error')
            return redirect(url_for('index'))
        time.sleep(5)
    return 'success'

@app.route('/IgG1_results/<id>')
def IgG1_results(id):
    return render_template('IgG1_result.html', IgG1_result=lines)


	#----------------------IgG3-------------------------

@app.route('/IgG3_predict',methods=['POST'])
def IgG3_predict():
    global filename
    global IgG3_email
    IgG3_email=request.form['IgG3_email'].replace(" ","")
    if request.form['IgG3_sequences'].replace(" ","")!="":
        filename=hashlib.md5(time.asctime()).hexdigest()
        file_in=open(filename, 'w')
        file_in.write(request.form['IgG3_sequences'].replace(" >",">"))
        file_in.close()
    else:
        file = request.files['IgG3_file']
        filename = secure_filename(file.filename)+"_"+hashlib.md5(time.asctime()).hexdigest()
        file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))

    task=run_IgG3.delay(filename, IgG3_email)
    id=task.id
    return redirect(url_for('IgG3_progress', id=id))



@app.route('/IgG3_progress/<id>')
def IgG3_progress(id):
    return render_template('IgG3_progress.html', id=id)
      
@app.route('/IgG3_status/<id>')
def IgG3_status(id):
    """ Wait until the job is finished and report success."""
    global lines
    result = run_IgG3.AsyncResult(id, app=celery)
    lines=result.get()
    while not result.ready():
        if result.failed():
            flash(u"Error running SchistoProt. please check input file", 'error')
            return redirect(url_for('index'))
        time.sleep(5)
    return 'success'

@app.route('/IgG3_results/<id>')
def IgG3_results(id):
    return render_template('IgG3_result.html', IgG3_result=lines)
	
	
#----------------------IgG4-------------------------

@app.route('/IgG4_predict',methods=['POST'])
def IgG4_predict():
    global filename
    global IgG4_email
    IgG4_email=request.form['IgG4_email'].replace(" ","")
    if request.form['IgG4_sequences'].replace(" ","")!="":
        filename=hashlib.md5(time.asctime()).hexdigest()
        file_in=open(filename, 'w')
        file_in.write(request.form['IgG4_sequences'].replace(" >",">"))
        file_in.close()
    else:
        file = request.files['IgG4_file']
        filename = secure_filename(file.filename)+"_"+hashlib.md5(time.asctime()).hexdigest()
        file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))

    task=run_IgG4.delay(filename, IgG4_email)
    id=task.id
    return redirect(url_for('IgG4_progress', id=id))



@app.route('/IgG4_progress/<id>')
def IgG4_progress(id):
    return render_template('IgG4_progress.html', id=id)
      
@app.route('/IgG4_status/<id>')
def IgG4_status(id):
    """ Wait until the job is finished and report success."""
    global lines
    result = run_IgG4.AsyncResult(id, app=celery)
    lines=result.get()
    while not result.ready():
        if result.failed():
            flash(u"Error running SchistoProt. please check input file", 'error')
            return redirect(url_for('index'))
        time.sleep(5)
    return 'success'

@app.route('/IgG4_results/<id>')
def IgG4_results(id):
    return render_template('IgG4_result.html', IgG4_result=lines)
	
	

app.wsgi_app = ProxyFix(app.wsgi_app)
if __name__ == "__main__":
    app.run(host='0.0.0.0', port=80, debug=True)# threaded=True)
