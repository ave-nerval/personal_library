# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

from dash import Dash, dcc, html, Input, Output, dash_table
import dash

import plotly.express as px
import pandas as pd

app = Dash(__name__)

# assume you have a "long-form" data frame
# see https://plotly.com/python/px-arguments/ for more options
df = pd.read_excel ('path')

meseci = ["Januar", "Februar", "Marec", "April", "Maj", "Junij", "Julij", "Avgust", "September", "Oktober", "November", "December"]

# popravi kasneje (dodaj še leto)
meseci_izbira = [x for x in meseci if x in df["Mesec"].unique()]
meseci_dropdown = []
for mesec in meseci_izbira:
    meseci_dropdown.append({'label': mesec, 'value': mesec})

df_vsota = df.groupby("Mesec")["Cena"].sum().round(2)
df_vsota2 = pd.DataFrame(data = {"Mesec": df_vsota.index.tolist(), "Vsota": df_vsota})


fig_vsota = px.scatter(x=df_vsota.index, y=df_vsota)
fig_vsota.update_layout(
    xaxis_title="Mesec",
    yaxis_title="Vsota"
    )


app.layout = html.Div(children=[
    html.H1(children='BAD HOBBITS'),

    html.Div(children='''
        Eva & Jernej
    '''),



    html.Div([
        dcc.Dropdown(
        options = meseci_dropdown, value = meseci_dropdown[0]["value"], multi = True, id='mesec-dropdown',
        style={'width': '30%', 'display': 'inline-block'}),


        dcc.RadioItems(['Vsota', 'Posebej'], 'Vsota', id = "tip-grafa", style={'display': 'inline-block'}),

    
    ],
    ),


    html.Div(children=[
        dcc.Graph(
            id='mesec-graf',
            style={'display': 'inline-block'}
        ),
        dcc.Graph(
            id='kategorija-graf',
            style={'display': 'inline-block'}
        )
    ]),


    html.Div(children=[
        dcc.Graph(
            id='vsota-graf',
            figure = fig_vsota,
            style={'display': 'inline-block'}
        ),



    ]),

    html.H3(children='Tabela (vsota)'),

    html.Div([
        dcc.Dropdown(
        options = ["Mesec", "Trgovina", "Kategorija"], value = "Mesec", multi = True, id='kategorija-dropdown'),
    
    ],
    style={'width': '90%', 'display': 'inline-block'}
    ),


    dash_table.DataTable(
            id = "table",
            style_data={
                'whiteSpace': 'normal',
                'height': 'auto',
            },
            fill_width=False,
            sort_action='native',
            filter_action="native",
            page_action='native',

        ),




])


@app.callback(
    [Output('mesec-graf', 'figure'),
    Output('kategorija-graf', 'figure')],
    [Input('mesec-dropdown', 'value'),
    Input('tip-grafa', 'value')]
)
def update_graph_trgovina(mesec, graf):
    if mesec:
        if isinstance(mesec, list):
            tmp = df.loc[df["Mesec"].isin(mesec)].copy()
            mesec_ime = str(mesec)[2:-2].replace("'", "")

        else:
            tmp = df.loc[df["Mesec"] == mesec].copy()
            mesec_ime = mesec

        if graf == "Vsota":
            fig = px.bar(tmp, x='Trgovina', y='Cena', title = "Mesec: " + mesec_ime, color = "Mesec", hover_data=["Opomba"])
            fig.update_layout(xaxis={'categoryorder':'total descending'})

            fig_kategorija = px.bar(tmp, x='Kategorija', y='Cena', color = "Mesec")
            fig_kategorija.update_layout(xaxis={'categoryorder':'total descending'})

        elif graf == "Posebej":
            fig = px.bar(tmp, x='Trgovina', y='Cena', title = "Mesec: " + mesec_ime, color = "Mesec", barmode="group")
            fig.update_layout(xaxis={'categoryorder':'total descending'})

            fig_kategorija = px.bar(tmp, x='Kategorija', y='Cena', color = "Mesec", barmode="group")
            fig_kategorija.update_layout(xaxis={'categoryorder':'total descending'})


        return fig, fig_kategorija
    
    else:
        return dash.no_update, dash.no_update
    


@app.callback(
    [Output("table", "data"), Output('table', 'columns')],
    Input('kategorija-dropdown', 'value')
)
def update_tabela_vsota(kategorija):
    if kategorija:
        tmp = df.groupby(kategorija)["Cena"].sum().round(2)
        tmp = tmp.to_frame()
        tmp = tmp.reset_index()
        cols = [{"name": i, "id": i} for i in tmp.columns]
        # tmp = tmp.to_dict('records'), [{"name": i, "id": i} for i in tmp.columns]
        return tmp.to_dict('records'), cols
    else:
        return dash.no_update   




if __name__ == '__main__':
    app.run_server(debug=True)



