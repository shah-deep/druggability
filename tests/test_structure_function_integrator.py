import pytest # type: ignore
import os
import json
from pathlib import Path
import sys

# Add project root to path to allow module imports
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from modules.structure_function_integrator import StructureFunctionIntegrator, ProteinMPNNResult

@pytest.fixture
def integrator():
    """Fixture to create and cleanup StructureFunctionIntegrator instance"""
    integrator = StructureFunctionIntegrator()
    yield integrator
    integrator.cleanup()

@pytest.fixture
def pdb_file():
    """Fixture to provide the test PDB file path"""
    return "examples/protein_42.pdb"

@pytest.fixture
def output_file():
    """Fixture to provide the test output file path"""
    return "test_sfi_results.json"

@pytest.fixture(autouse=True)
def setup_teardown(output_file):
    """Fixture to handle setup and teardown for tests"""
    # Ensure we are in the project root for file paths to work
    # This is important for locating the example pdb file
    os.chdir(Path(__file__).resolve().parents[1])
    yield
    # Cleanup after each test
    if os.path.exists(output_file):
        os.remove(output_file)

def test_analyze_structure_mock_mode(integrator, pdb_file, output_file):
    """Test the analyze_structure method in mock mode"""
    # Since ProteinMPNN is not installed, this will run the mock analysis
    result = integrator.analyze_structure(pdb_file)

    assert isinstance(result, ProteinMPNNResult)
    assert isinstance(result.binding_site_score, float)
    assert 0.01 <= result.binding_site_score <= 1.0
    assert "ProteinMPNN" in result.model_used

    # Test exporting the results
    integrator.export_results(result, output_file)
    assert os.path.exists(output_file)
    
    with open(output_file, 'r') as f:
        data = json.load(f)
    assert data['binding_site_score'] == result.binding_site_score
    assert "model_used" in data

def test_protein_mpnn_result_creation():
    """Test ProteinMPNNResult creation and attributes"""
    result = ProteinMPNNResult(
        binding_site_score=0.75,
        model_used="ProteinMPNN_v1.0",
        sequence_recovery=0.85,
        structural_plausibility=0.92,
        druggability_metrics={"confidence": 0.88, "stability": 0.91},
        designed_sequences=["MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"],
        confidence_scores=[0.92, 0.88, 0.95],
        warnings=["Low confidence in binding site prediction"]
    )
    
    assert result.binding_site_score == 0.75
    assert result.sequence_recovery == 0.85
    assert result.model_used == "ProteinMPNN_v1.0"
    assert result.structural_plausibility == 0.92
    assert len(result.designed_sequences) == 1
    assert len(result.confidence_scores) == 3
    assert result.druggability_metrics["confidence"] == 0.88
    assert result.druggability_metrics["stability"] == 0.91
    assert len(result.warnings) == 1

def test_integrator_initialization(integrator):
    """Test that the integrator initializes correctly"""
    assert integrator is not None
    assert hasattr(integrator, 'analyze_structure')
    assert hasattr(integrator, 'export_results')
    assert hasattr(integrator, 'cleanup') 