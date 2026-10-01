#!/bin/bash
# Launch the TRADES Config Generator GUI
cd "$(dirname "$0")"
streamlit run app.py "$@"
