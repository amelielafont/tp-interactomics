# tp-interactomics


## Contexte biologique

Nous allons reproduire certaines analyses rapportées dans une étude du réseau d'interaction des protéines du [virus d'Epstein-Bar](https://en.wikipedia.org/wiki/Epstein%E2%80%93Barr_virus) avec certaines protéines de son hôte principal l'homme ([Calderwood et al.](https://www.pnas.org/content/104/18/7606)). Le cible cellulaire principale du virus EBV est lymphocyte B humain.  

## Mise en place
Vous commencerez par "fourcher" et cloner ce repository.
Seul Jupyter avec la librarire [networkx](https://networkx.org/) est requis.
Les libraries Pandas et request peuvent également être utilisées.
A l'exception des encarts du présent `README` prévus à cet effet, il vous est conseillé de réaliser le TP dans un notebook Jupyter que vous joindrez à ce répository git en fin de séance.

### Données

#### Protéomiques

Les fiches UNIPROT des protéines étudiées dans la publication sont disponibles dans fichiers XML suivants:

- `data/Calderwood_EBV_proteome.xml`
- `data/Calderwood_Human_proteome.xml`

#### Interactomiques

Nous allons devoir recupérer les données d'interactions étudiées dans la publication grâce au protocole [PSICQUIC](https://psicquic.github.io/PsicquicSpec_1_4_Rest.html).

Ce protocole permet l'accès à distance à de nombreuses bases de données d'interactions protéine-protéine. Les interactions présentes dans ces bases de données sont obtenues par curation minutieuse de la littérature scientifique. Les interactions sont uniquement binaires (2 protéines) et toujours associées à la publication d'origine.
D'autres informations peuvent également être rapportées.
Le format standard de PSICQUIC, tabulé, est nommé [MITAB](https://psicquic.github.io/PSIMITAB.html).
Les requêtes sont formulées aux bases de données à l'aide d'un protocole REST obéissant à la syntaxe [MIQL](https://psicquic.github.io/MiqlReference.html)

Pour ce TP nous utiliserons la bases de données [Intact](https://www.ebi.ac.uk/intact/), dont le service PSICQUIC est accessible à l'URL: `http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search`.
Par exemple, les 100 premières interactions protéine-protéine humaines disponibles dans Intact sont accessibles via l'URL: `http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/species:human?firstResult=0&maxResults=100`

##### Quelles sont les significations des champs suivants du format MITAB 2.8?

Numero de champ | Signification Biologique|
 --- | --- 
1 | Unique identifier for interactor A
2 | Unique identifier for interactor B
3 | Alternative identifier for interactor A
4 | Alternative identifier for interactor B
5 | Aliases for A
6 | Aliases for B
7 | Interaction detection methods
8 | First author
9 | Identifier of the publication

##### Utiliser le PMID(17446270) de la publication pour récuperer les lignes MITAB des interactions rapportées dans l'étude.
Une librairie pratique pour manipuler des requêtes HTTP est [requests](https://requests.readthedocs.io/en/master/), eg:

```python
import requests
url = "https://mmsb.cnrs.fr/equipe/mobi/"

try:
    httpReq = requests.get(url, proxies=None)
except NameError:
    httpReq = requests.get(url)
ans = httpReq.text
```
Out: 
```
uniprotkb:Q82235	uniprotkb:Q9DBG9	intact:EBI-1162486	intact:EBI-1161647	psi-mi:q82235_9dela(display_long)|uniprotkb:tax(gene name)|psi-mi:tax(display_short)	psi-mi:tx1b3_mouse(display_long)|uniprotkb:Tax interaction protein 1(gene name synonym)|uniprotkb:Tax1bp3(gene name)|psi-mi:Tax1bp3(display_short)	psi-mi:"MI:0096"(pull down)	Kanamori et al. (2003)	pubmed:12874278	taxid:11908(humt-)|taxid:11908(Human T-lymphotropic virus 1)	taxid:10090(mouse)|taxid:10090(Mus musculus)	psi-mi:"MI:0407"(direct interaction)	psi-mi:"MI:0469"(IntAct)	intact:EBI-1162495	intact-miscore:0.44
uniprotkb:P04578	uniprotkb:B0FAM1	intact:EBI-6163496|uniprotkb:O09779	intact:EBI-6163546	psi-mi:env_hv1h2(display_long)|uniprotkb:env(gene name)|psi-mi:env(display_short)|uniprotkb:Env polyprotein(gene name synonym)	psi-mi:b0fam1_9hiv1(display_long)|uniprotkb:env(gene name)|psi-mi:env(display_short)	psi-mi:"MI:0007"(anti tag coimmunoprecipitation)	Jager et al. (2012)	psi-mi:"MI:0007"|imex:IM-17346|pubmed:22190034	taxid:11706(hv1h2)|taxid:11706("Human immunodeficiency virus type 1 group M subtype B (isolate HXB2) (HIV-1)")	taxid:11676(9hiv1)|taxid:11676(Human immunodeficiency virus 1)	psi-mi:"MI:0914"(association)	psi-mi:"MI:0469"(IntAct)	intact:EBI-6174316|imex:IM-17346-4	intact-miscore:0.35
.
.
.
uniprotkb:Q9Q2G4	uniprotkb:Q9UGP4	intact:EBI-6248094	intact:EBI-2652871|uniprotkb:Q17RQ1|uniprotkb:Q9BQQ9|uniprotkb:Q9NQ47	psi-mi:q9q2g4_9hiv1(display_long)|uniprotkb:ORF(gene name)|psi-mi:ORF(display_short)|uniprotkb:Retropepsin(gene name synonym)	psi-mi:limd1_human(display_long)|uniprotkb:LIMD1(gene name)|psi-mi:LIMD1(display_short)	psi-mi:"MI:0096"(pull down)	Jager et al. (2012)	psi-mi:"MI:0007"|imex:IM-17346|pubmed:22190034	taxid:11676(9hiv1)|taxid:11676(Human immunodeficiency virus 1)	taxid:9606(human)|taxid:9606(Homo sapiens)	psi-mi:"MI:0914"(association)	psi-mi:"MI:0469"(IntAct)	intact:EBI-6177860|imex:IM-17346-61	intact-miscore:0.56

```

##### Quelles techniques experimentales mesurent les interactions rapportées dans cette publication?

```

```

#### Extraction des deux sous-jeux d'interactions

Vous disposez d'un code écrit par votre collègue pour extraire à partir d'un chaîne de caractères au format mitab les lignes concertant:

- Les interactions EBV-EBV
- Les interactions EBV-Humaine

Malheuresement, votre collègue a du partir en urgence. Vous devez comprendre et compléter son code pour mener à bien l'extraction des deux jeux d'interactions.

```python
import re

def mitabReader(httpText):
    for line in ans.split("\n"):
        _ = line.split("\t")
        if len(_) > 1 and _[0].startswith("uniprotkb:")\
                      and _[1].startswith("uniprotkb:"):
            yield [ _[0].replace("uniprotkb:", ""),\
                    _[1].replace("uniprotkb:", "") ]\
                  + _[2:]
                
            
def isMitab_EBV_EBV(mitabArray):
    reEBV   = "taxid:(1037[6-7]|82830)"
    if re.search(reEBV, mitabArray[9]) and re.search(reEBV, mitabArray[10]):
        return True
    return False

def isMitab_Human_EBV(mitabLine):
    reHuman = 'taxid:9606'
    if re.search(reHuman, mitabArray[9]) or re.search(reHuman, mitabArray[10]):
        return True
    return False


EBV_EBV_mitab   = []
EBV_Human_mitab = []
total = 0
for mitabArray in mitabReader(ans):
    total += 1
    if isMitab_EBV_EBV(mitabArray):
        EBV_EBV_mitab.append(mitabArray)
    elif isMitab_Human_EBV(mitabArray):
        EBV_Human_mitab.append(mitabArray)
    else : 
        raise ValueError("Je ne connais pas cette espece ==> ", mitabArray[9:11])

print(f"Nombre total d'interactions {total}, EBV-EBV {len(EBV_EBV_mitab)}")
```

##### Que fait la fonction `mitabReader` ?
```
```

##### Après avoir réparé ce code veuillez
- Extraire les lignes MITAB impliquant uniquement des protéines d'EBV, quel est leur nombre ?
- Extraire les lignes MITAB impliquant des protéines humaines et des protéines d'EBV, quel est leur nombre ?
```
Nombre total d'interactions 230, EBV-EBV 59, EBV-Human 171
```

##### Combien de protéines humaines et virales sont respectivement dans les jeux d'interactions EBV-Human et EBV-EBV ?

```python
EBV_protein = set()
for data in EBV_EBV_mitab:
    EBV_protein.add(data[0])
    EBV_protein.add(data[1])
print(f"{len(EBV_protein)} EBV proteins")

Human_protein = set()
for data in EBV_Human_mitab:
    Human_protein.add(data[0])
    Human_protein.add(data[1])

Human_protein = Human_protein - EBV_protein
print(f"{len(Human_protein)} Human proteins")
```
Out: 
```
48 EBV proteins
129 Human proteins
```

###### Pour la suite du travail assurez-vous d'avoir les deux jeux de données MITAB suivants

- MITAB EBV/EBV
- MITAB EBV/Human

Chacun ne contetant que des interactants référés par un numéro d'accession UNIPROT

### Construction du réseau d'interactions EBV/EBV

A l'aide des données MITAB et de la librarie [networkx](https://networkx.github.io/documentation/latest/
), représentez graphiquement un réseau dans lequel:

- les noeuds sont des identifiants UNIPROT

- les arêtes relient deux protéines en interaction

```python
%matplotlib inline
import matplotlib.pyplot as plt

import networkx as nx
G = nx.Graph()

plt.figure(figsize=(10,8))

for data in EBV_EBV_mitab:
    p1 = data[0]
    p2 = data[1]
    G.add_edge(p1, p2)

nx.draw(G, with_labels=True, font_weight='bold', node_size=500)
#plt.show()
#or
plt.savefig("ebv_ebv_network_uniprot.png")
```

![Graphique](MADP_TP2/ebv_ebv_network_uniprot.png)

##### Décrivez brièvement ce réseau

```
taille de la plus grande composante, discussion avec le nom des gènes.
```

Les noms de gènes sont parfois plus parlants que des accesseurs UNIPROT. A l'aide du fichier `./data/Calderwood_EBV_proteome.xml`  créez une table de conversion entre accesseur UNIPROT et nom de gène.

Pour vous aidez dans cette tâche, vous disposez du "parser" XML suivant qui étant donné un numéro d'accession UNIPROT et un document XML retourne un dictionaire d'informations concernant cette protéine. Le code permettant d'extraire l'information du nom de gène à été supprimée lors d'un `git push` malencontreux, à vous de le completer avant utilisation.

```python
from xml.etree.ElementTree import parse, dump, fromstring, register_namespace, ElementTree

# Utility functions
# Extracting All go terms relative to provided UNIPROT accessor
def goTerms(xmlEntry):
    ns = '{http://uniprot.org/uniprot}'
    goTerms = xmlEntry.findall(ns +'dbReference[@type="GO"]')
    goTermList = []
    for goT in goTerms:
        gID   = goT.attrib['id']
        gName = goT.find(ns +'property[@type="term"]').attrib['value']
        goTermList.append({"name" : gName, "ID" : gID})
    return goTermList

# Return information about provided UNIPROT accessor as python dictionary
def proteinDict(uniprotID, root):
    ns   = '{http://uniprot.org/uniprot}'

    data = { "accession" : uniprotID,
             "geneName" : None,
             "name" : None,
             "GOterms" : None
           }

    for entry in root.findall(ns+'entry'):
        accessions = entry.findall(ns+"accession")
        for acc in accessions:
            if acc.text == uniprotID: # entry is the node matching provided UNIPROT accessor
                e = entry.find(f"{ns}protein/{ns}recommendedName/{ns}fullName")
                if not e is None:
                    data["name"] = e.text
                e = "OUPSS##!!!"
                if not e is None:
                    data["geneName"] = e.text

                data["GOterms"] = goTerms(entry)
                return data
    raise ValueError(f"{uniprotID} nor found in XML document")
```

```python
# Test
tree = parse('./data/Calderwood_Human_proteome.xml')
root = tree.getroot()
proteinDict("Q53Y88", root)
```

Vous pouvez desormais dessiner le réseau dans lequel:

- les arêtes relient deux protéines en interaction
- les noeuds sont les noms des gènes correspondant aux protéines.

Correction:
```python
tree = parse('../data/Calderwood_EBV_proteome.xml')
root = tree.getroot()

EBV_node_label = {
    #uniprotID : GeneName // uniprotID
}

for d in EBV_EBV_mitab:
    for uniprotID in d[:2]:
        if not uniprotID in EBV_node_label:
            protData = proteinDict(uniprotID, root)
            EBV_node_label[uniprotID] = protData['geneName'] if protData['geneName'] else uniprotID
EBV_node_label

%matplotlib inline
import matplotlib.pyplot as plt

import networkx as nx
G = nx.Graph()

plt.figure(figsize=(10,8))

for data in EBV_EBV_mitab:
    p1 = data[0]
    p2 = data[1]
    G.add_edge(p1, p2)

G = nx.relabel_nodes(G, EBV_node_label)
nx.draw(G, with_labels=True, font_weight='bold', node_size=500)
#plt.show()
#or
plt.savefig("ebv_ebv_network_gene.png")
```
![Graphique](MADP_TP2/ebv_ebv_network_gene.png)

### Caractérisation des cibles protéiques du virus

Pour perturber et détourner à son profit le fonctionnement du lymphocyte, le matériel protéique du virus va interagir préferentiellement avec certains processus cellulaires.

En suivant, la méthode précédente dessiner le réseau dans lequel:

- les arêtes relient une protéine humaine et une protéine virale en interaction.
- les noeuds sont les noms des gènes correspondant aux protéines.
- les noeuds des protéines virales et humaines sont dessinés différement.

Vous pouvez jouer sur la taille de la figure et la constante de ressort *k* du rendu [spring_layout](https://networkx.org/documentation/stable/reference/generated/networkx.drawing.layout.spring_layout.html) pour augmenter la lisibilité du graphique

Correction :
```python
plt.figure(figsize=(24,18))
G =nx.Graph()

for data in EBV_Human_mitab:
    G.add_edge(data[0], data[1])

#get human and ebv node list
humanNodeName = list(Human_node_label.keys())
ebvNodeName = list( set(list(G.nodes())) - set(humanNodeName))

#compute layout then draw on it
pos = nx.spring_layout(G, k=0.1)

#node fonts, shapes and colors
nx.draw_networkx_nodes(G, pos, ebvNodeName, node_color='red', node_shape='s')
nx.draw_networkx_nodes(G, pos, humanNodeName, node_color='blue', node_shape='o')
nx.draw_networkx_edges(G, pos)
nx.draw_networkx_labels(G, pos, labels=Human_node_label, font_size = 12)

#plt.show()
#or
plt.savefig("ebv_human_network_gene.png")
```

![Graphique](MADP_TP2/ebv_human_network_gene.png)








## Analyse des interactions entre protéines humaines avec STRING

Maintenant que vous avez pu visualiser le réseau des interactions physiques virus/homme, vous allez utiliser la ressource STRING pour visualiser les interactions entre les protéines humaines ciblées par le virus.

Comme décrit dans la publication, 173 interactions virus/homme ont été observées, impliquant 112 protéines humaines.


### Representation des interactions avec STRING

Aller sur le site de [STRING](string-db.org), et afficher les interactions entre les 112 protéines humaies ciblées par le virus:

1. Aller dans *Search*
2. puis *Multiple proteins*
3. charger le fichier human_proteins_unique.txt, en spécifiant l'organisme *Homo Sapiens*,
4. valider les protéines identifiées par STRING (en cliquant sur *Continue*).

Par défaut, STRING affiche les interactions avec tous les canaux, a partir d'un indice de confiance 'moyen', c'est a dire supérieur a 0.4.
Nous allons conserver ces paramètres par défaut.

## Simplification du réseau

Nous allons essayer d'exploiter ce réseau pour mettre en évidence les voies fonctionnelles ciblées par le virus.

Afin de simplifier l'analyse, éliminer les protéines non connectées au reste du réseau:
Aller dans Settings, Advanced Settings, network display options: hide disconnected nodes in the network.

Que pouvez-vous dire de la structure de ce réseau ?
Combien de 'modules' (= clusters) semble-t-il y avoir ?

Confirmez votre interprétation en vous aidant de l'outil de clustering de STRING.
Uploadez une image du réseau coloré qui vous semble le mieux représenter sa structure interne.

![RESEAU_STRING](ebv_ebv_network_gene.png)

## Identification des voies fonctionnelles

Pour les deux clusters principaux, naviguez dans les annotations des protéines, et essayez de proposer des hypothèse quand au mécanisme d'action de EBV, en expliquant votre démarche.

(Toutes les hypothèses sont intéressantes si vous expliquez votre façon de faire).

```
```

Comment procéderiez-vous si vous deviez faire cette analyse de manière plus rigoureuse ?

```

```







## Questions Optionelles

Les noeuds du réseau d'interaction peut aussi être divisés en deux partitions, humaine et virale. Vous disposez, ci-dessous d'un exemple de rendu graphique "bipartite" sur deux selections arbitraires de noeuds. Essayez de l'adapter au problème de représentation graphique du réseau d'interaction EBV/Humain

```python
import networkx as nx
plt.figure()
plt.axis('off')

G = nx.Graph()
G.add_edge("A", "C")
G.add_edge("B", "C")
X, Y = nx.bipartite.sets(G,  ["A", "B"])
pos = nx.bipartite_layout(G, X)
nx.draw_networkx_nodes(G, pos, node_color="blue", node_shape="o",nodelist=["A", "B"] )
nx.draw_networkx_nodes(G, pos, node_color="red", node_shape="o",nodelist=["C"] )
nx.draw_networkx_labels(G,pos,font_weight=800,font_color='black')
nx.draw_networkx_edges(G, pos, width=0.5)
```

##### Propriétés topologiques des protéines humaines ciblées par le virus

Les paires connues de protéines humaines en interaction sont disponibles dans le fichier `data/human_interactome_pair_uniprot.tab.gz`.

* Construisez l'objet `networkx.Graph` stockant ces interactions.
* Dessinez l'histogramme multiple des distributions des degrés dans cet objet réseau pour les deux populations de protéines suivantes:
   - Protéines humaines ciblées par EBV
   - Protéines humaines totales

![Graphique](human_protein_degree.png)

Qu'observez vous ?
```
```

##### Identification des processus biologiques ciblés par le virus

Nous allons observés les termes GO présents dans les protéines humaines les plus ciblées par EBV

1. Classer les noeuds, protéines humaines, par degré décroissant.
2. Lister les termes GO des protéines les plus connectées

- Quels termes reviennent fréquemment?

```

```

- Connaissant l'organisation de l'ontologie GO, quelles suggestions pourriez vous faire afin d'augmenter la qualité de l'analyse?

```

```
