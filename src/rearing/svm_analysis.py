from sklearn import svm, preprocessing
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, classification_report

df = df_frametimeinfo.join(df_stimulus.set_index('stim_idx'),
                           on='peakresp_idx', how='inner')

y = df['stimname'].map({'1-6ol': -1, 'EP': 1}).dropna()
Y = df['stimname'].loc[y.index].to_numpy()
X = dfof[:, y.index].T
#%%
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=.25, random_state=42)

svc = svm.SVC(kernel='linear')
svc.fit(X_train, y_train)
y_pred = svc.predict(X_test)

# accuracy
accuracy = accuracy_score(y_train, svc.predict(X_train))
print(f"\naccuracy: {accuracy}")

# classification report
print(classification_report(y_test, y_pred))
#%%
# svm projection distance
svm_prj = svc.decision_function(dfof.T)

# plot projection distance to hyperplane
fig = px.line(x=df_frametimeinfo.frame_times, y=svm_prj)
fig.show()

fig.add_vrect(x0="2018-09-24", x1="2018-12-18", col=1,
              annotation_text="decline", annotation_position="top left",
              fillcolor="green", opacity=0.25, line_width=0)

# %%
classifiers = {
    'Linear SVC': svm.SVC(kernel='linear', C=C, probability=True,
                          random_state=0),

}

n_classifiers = len(classifiers)

# %%
svc = svm.SVC(kernel='linear')
svc.fit(scaler.transform(X_train), y_train)

dfof_prj = svc.decision_function(scaler.transform(dfof.transpose()))

svc = svm.SVC(kernel='linear')
svc.fit(X_train, y_train)
# %%
y_pred = svc.predict(X_train)
accuracy = accuracy_score(y, y_pred)
print("Accuracy (train) for %s: %0.1f%% " % (name, accuracy * 100))
# %%
svc = svm.SVC(kernel='linear')
svc.fit(scaler.transform(X_train), y_train)
svm_prj = svc.decision_function(scaler.transform(dfof.transpose()))

fig = px.line(x=df_frametimeinfo.frame_times, y=dfof_prj);
fig.show()
#%%
fig.add_vrect(x0="2018-09-24", x1="2018-12-18", col=1,
              annotation_text="decline", annotation_position="top left",
              fillcolor="green", opacity=0.25, line_width=0)
