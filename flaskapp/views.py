from flaskinsight import app
from flask import render_template, redirect, request
import pandas as pd
import os, sys
import pickle
import sklearn.linear_model as sk

filepath = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(filepath, "models", "Targets.pickle"), "rb") as f:
    targets = pickle.load(f)
with open(os.path.join(filepath, "models", "SideEffects.pickle"), "rb") as f:
    sideeffects = pd.DataFrame(pickle.load(f))
def conf(preds,thr = 0.5):
    clf = preds.copy(deep=True)
    clf[clf >=thr] = 1
    clf[clf <= thr] = 0
    clf = clf.sum(axis = 1)
    return clf
def risks(preds,thr = 0.5):
    clf = preds.copy(deep = True)
    clf[clf <= thr] = 0.
    clf = clf.sum(axis = 1)
    return clf

@app.route('/')
def index():
    return render_template("starter-template.html")

@app.route('/about')
def about():
    return render_template("about.html")

@app.route('/contact')
def contact():
    return render_template("contact.html")

@app.route('/', methods = ["POST"])
def result():
    #Open features
    with open(os.path.join(filepath, "models", "Targets.pickle"), "rb") as f:
        targets = pickle.load(f)
    targets = pd.DataFrame(targets)

    #read in the inputs
    molecule1 = request.form['molecule1']
    if molecule1 == "":
        molecule1 = "Drug1"
    target1 = request.form['target1']

    molecule2 = request.form['molecule2']
    if molecule2 == "":
        molecule2 = "Drug2"
    target2 = request.form['target2']

    molecule3 = request.form['molecule3']
    if molecule3 == "":
        molecule3 = "Drug3"
    target3 = request.form['target3']

    select_n = int(request.form['select_n'])
    if select_n == "":
        select_n = 1

    #Make feature matrix

    drugs = [molecule1, molecule2, molecule3]
    inputs = [target1, target2, target3]

    X_test = pd.DataFrame(index = drugs, columns = targets)
    X_input = []
    index = 0

    for inp in inputs:
        targ = inp.split("\r\n")
        targ = [m.lower() for m in targ]

        X_input = pd.DataFrame([int(t in targ) for t in targets.iloc[:,0].str.lower().values])
        X_test.loc[drugs[index],:] = pd.DataFrame.transpose(X_input).values
        index += 1


    #Given the models we made above, predict output for user-provided features
    #Use model to predict each side effect


    #Initialize predictions
    cols = ['Model1','Model2','Model3','Model4','Model5']

    preds_nausea = pd.DataFrame(index = drugs, columns = cols)
    preds_fever = pd.DataFrame(index = drugs, columns = cols)
    preds_dyspnoea = pd.DataFrame(index = drugs, columns = cols)
    preds_rash = pd.DataFrame(index = drugs, columns = cols)
    preds_vomiting = pd.DataFrame(index =  drugs, columns = cols)
    m_number = 0

    #Loop through models
    for m_number in range(0,len(cols)):
        #Pre-process feature input.
        with open(os.path.join(filepath, "models", 'poly_%s.pickle' % m_number), 'rb') as f:
            poly = pickle.load(f)
        with open(os.path.join(filepath, "models", 'sfm_%s.pickle' % m_number), 'rb') as f:
            sfm = pickle.load(f)
        with open(os.path.join(filepath, "models", 'Targets.pickle'), 'rb') as f:
            targets = pickle.load(f)
        X_test_new = pd.DataFrame(poly.fit_transform(X_test))
        X_test_new.columns = poly.get_feature_names(targets)
        X_test_transform = pd.DataFrame(sfm.transform(X_test_new))
        X_test_transform.columns = X_test_new.columns[sfm.get_support()]

        #Open model and predict each side effects
        with open(os.path.join(filepath, "models",'LogisticRegression_%s.pickle' % m_number), 'rb') as f:
            model = pickle.load(f)

        preds_nausea.iloc[:,m_number] = pd.DataFrame(model.predict_proba(X_test_transform)).iloc[:,0].values
        preds_fever.iloc[:,m_number] = pd.DataFrame(model.predict_proba(X_test_transform)).iloc[:,1].values
        preds_dyspnoea.iloc[:,m_number] = pd.DataFrame(model.predict_proba(X_test_transform)).iloc[:,2].values
        preds_rash.iloc[:,m_number] = pd.DataFrame(model.predict_proba(X_test_transform)).iloc[:,3].values
        preds_vomiting.iloc[:,m_number] = pd.DataFrame(model.predict_proba(X_test_transform)).iloc[:,4].values
        m_number = m_number + 1

    #Combine model predictions into output prediction matrix.
    #Initialize matrix
    y=[]
    y = pd.DataFrame(index = X_test.index, columns = sideeffects.iloc[:,0])
    preds = [preds_nausea, preds_fever, preds_dyspnoea, preds_rash, preds_vomiting]
    index = 0

    for prd in preds:
        y.iloc[:,index] = prd.median(axis = 1)
        index += 1

    #Determine best performers
    y_clf = y.copy(deep = True)
    y_clf.loc[:,'clf'] = pd.DataFrame(conf(y), index = y.index)
    y_acc = y_clf.sort_values(by = 'clf', ascending = True).iloc[0:select_n,:]
    y_unacc = y_clf.sort_values(by = 'clf', ascending = True).iloc[select_n:,:]
    if min(y_unacc.loc[:,'clf']) < max(y_acc.loc[:,'clf']):
        y_suggested = y_acc
    else:
        y_proba = y.copy(deep = True)
        y_proba.loc[:,'total'] = pd.DataFrame(risks(y))
        y_acc = y_proba.sort_values(by = 'total', ascending = True).iloc[0:select_n,:]
        y_unacc = y_proba.sort_values(by = 'total', ascending = True).iloc[select_n:,:]
        y_suggested = y_acc



    #Convert predictions to strings
    drugcandidates = "We recommend proceeding with " + ', '.join(y_suggested.index)
    pred_string = []
    sentencestructure = pd.DataFrame(index = y_suggested.index)
    for d in y_suggested.index:
        pred_string = []
        for s in sideeffects.iloc[:,0]:
            if y_suggested.loc[d,s] > 0.5:
                pred_string += [s]

        if not len(pred_string):
            sentencestructure.loc[d, 'String'] = "No adverse reactions anticipated for " + d
        elif len(pred_string) > 1:
            sentencestructure.loc[d,'String'] = d + " may cause " + ', '.join(neg_pred_string[:-1]) +" and " +neg_pred_string[-1]
        else:
            sentencestructure.loc[d,'String'] = d + " may cause " + ', '.join(pred_string)
    prediction = '. '.join(sentencestructure.iloc[:,0])

    #Explain negative predictions
    negdrugcandidates = "We do not recommend proceeding with " + ', '.join(y_unacc.index)
    neg_pred_string = []
    negsentencestructure = pd.DataFrame(index = y_unacc.index)
    for d in y_unacc.index:
        neg_pred_string = []
        for s in sideeffects.iloc[:,0]:
            if y_unacc.loc[d,s] > 0.5:
                neg_pred_string += [s]

        if not len(neg_pred_string):
            negsentencestructure.loc[d, 'String'] = "No adverse reactions anticipated for " + d
        elif len(neg_pred_string) > 1:
            negsentencestructure.loc[d,'String'] = d + " has a higher likelihood of causing " + ', '.join(neg_pred_string[:-1]) +" and " +neg_pred_string[-1]
        else:
            negsentencestructure.loc[d,'String'] = d + " has a higher likelihood of causing " + ', '.join(neg_pred_string)
    negprediction = '. '.join(negsentencestructure.iloc[:,0])

    return render_template("result-template.html", drugcandidates = drugcandidates, prediction = prediction, sentencestructure = sentencestructure, \
    negprediction=negprediction, negdrugcandidates = negdrugcandidates, negsentencestructure = negsentencestructure, table = (100*y).to_html(float_format = '{:.0f}%'.format, classes = ["table", "table-striped"]))
