if [ ! -d ".env" ]; then
    echo "Setting up virtual environment for file conversion and activating it..."
    python3 -m venv .env && source .env/bin/activate
    pip install h5py uproot tqdm
else
    echo "Virtual environment for file conversion already set, activating it..."
    source .env/bin/activate
fi
