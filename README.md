# Vascularisé Fractal Volume Analysis

Projet Python pour analyser un volume vasculaire 3D à partir d'un masque binaire, extraire le squelette vasculaire, calculer les rayons locaux et estimer une dimension fractale ainsi qu'un volume sanguin total extrapolé.

compatibilité PARFAITE AVEC LA HAS ET LES PROTOCOLES D’ESSAI ET DE TESTS CLINIQUES POUR LA LUTTE CONTRE LES MALADIES VIRALES TYPE SIDA OU DE L’HERPÈS ET l'HEPATITE GENERIQUE PAR EXEMPLE.

## Description

Le script `notre_Vd.py` contient une classe `FractalVascularVolume` qui permet de :

- convertir un volume binaire en masque booléen,
- calculer la transformée de distance (EDT),
- extraire le squelette topologique 3D,
- calculer les rayons locaux du réseau,
- ajuster une loi de puissance en log-log,
- estimer la dimension fractale locale,
- extrapoler le volume microvasculaire et le volume total effectif.

## Installation

1. Créer un environnement virtuel :

```bash
python3 -m venv .venv
source .venv/bin/activate
```

2. Installer les dépendances :

```bash
pip install -r requirements.txt
```

## Utilisation

Exécuter le script principal :

```bash
python notre_Vd.py
```

Le script génère un volume vasculaire factice pour tester le pipeline complet.

## Dépendances

- numpy
- scipy
- scikit-image

## Structure du projet

```text
Vd/
├── .gitignore
├── README.md
├── requirements.txt
├── notre_Vd.py
└── .venv/   # ignoré par Git
```

## Licence

Ce projet est fourni à titre éducatif et peut être adapté selon vos besoins.
