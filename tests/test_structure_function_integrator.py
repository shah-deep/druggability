import pytest # type: ignore
import os
import json
from pathlib import Path
import sys
from unittest.mock import patch, MagicMock
# Mock numpy for testing
class MockNumpy:
    def random(self):
        return MockNumpy()
    
    def __call__(self, *args, **kwargs):
        return [[0.1] * 21] * 100  # Mock probabilities

np = MockNumpy()

# Add project root to path to allow module imports
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from modules.structure_function_integrator import StructureFunctionIntegrator, ProteinMPNNResult

@pytest.fixture
def integrator():
    """Fixture to create and cleanup StructureFunctionIntegrator instance"""
    # Mock the model setup to avoid sys.exit(1)
    with patch('modules.structure_function_integrator.StructureFunctionIntegrator._setup_proteinmpnn_model') as mock_setup:
        mock_setup.return_value = "/mock/model/path.pt"
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
    """Test the analyze_structure method with mocked ProteinMPNN execution"""
    # Mock the ProteinMPNN execution to simulate successful run
    with patch.object(integrator, '_run_proteinmpnn') as mock_run:
        mock_run.return_value = {
            'sequences': ['MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG'],
            'scores': [0.92, 0.88, 0.95],
            'recovery': 0.85,
            'probabilities': [[[0.1] * 21] * 100] * 10  # Mock probabilities
        }
        
        # Mock binding site extraction
        with patch.object(integrator, '_extract_binding_sites') as mock_binding:
            mock_binding.return_value = [{'chain': 'A', 'center': [0, 0, 0], 'residue_range': [1, 10], 'volume': 100}]
            
            result = integrator.analyze_structure(pdb_file)

            assert isinstance(result, ProteinMPNNResult)
            assert isinstance(result.binding_site_score, float)
            assert 0.0 <= result.binding_site_score <= 1.0  # Allow full range
            assert "ProteinMPNN_v_48_020" in result.model_used  # Correct format
            assert result.sequence_recovery == 0.85
            assert len(result.designed_sequences) == 1
            assert len(result.confidence_scores) == 3

            # Test exporting the results
            integrator.export_results(result, output_file)
            assert os.path.exists(output_file)
            
            with open(output_file, 'r') as f:
                data = json.load(f)
            assert data['binding_site_score'] == result.binding_site_score
            assert "model_used" in data
            assert "timestamp" in data
            assert "config" in data

def test_protein_mpnn_result_creation():
    """Test ProteinMPNNResult creation and attributes"""
    result = ProteinMPNNResult(
        binding_site_score=0.75,
        model_used="ProteinMPNN_v_48_020",
        sequence_recovery=0.85,
        structural_plausibility=0.92,
        druggability_metrics={"confidence": 0.88, "stability": 0.91},
        designed_sequences=["MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"],
        confidence_scores=[0.92, 0.88, 0.95],
        warnings=["Low confidence in binding site prediction"]
    )
    
    assert result.binding_site_score == 0.75
    assert result.sequence_recovery == 0.85
    assert result.model_used == "ProteinMPNN_v_48_020"
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
    assert hasattr(integrator, 'config')
    assert hasattr(integrator, 'logger')

def test_pdb_validation(integrator):
    """Test PDB structure validation"""
    # Test with valid PDB file
    valid_pdb = "examples/protein_42.pdb"
    assert integrator._validate_pdb_structure(valid_pdb) == True
    
    # Test with non-existent file
    with pytest.raises(FileNotFoundError):
        integrator._validate_pdb_structure("nonexistent.pdb")
    
    # Test with empty file
    empty_file = "test_empty.pdb"
    try:
        with open(empty_file, 'w') as f:
            f.write("")
        assert integrator._validate_pdb_structure(empty_file) == False
    finally:
        if os.path.exists(empty_file):
            os.remove(empty_file)

def test_error_handling(integrator, pdb_file):
    """Test error handling when ProteinMPNN fails"""
    # Mock ProteinMPNN to raise an exception
    with patch.object(integrator, '_run_proteinmpnn') as mock_run:
        mock_run.side_effect = RuntimeError("ProteinMPNN failed")
        
        with pytest.raises(RuntimeError):
            integrator.analyze_structure(pdb_file)

def test_export_results_format(integrator, output_file):
    """Test that exported results have the correct format"""
    result = ProteinMPNNResult(
        binding_site_score=0.8,
        model_used="ProteinMPNN_v_48_020",
        sequence_recovery=0.9,
        structural_plausibility=0.85,
        druggability_metrics={"confidence": 0.9, "stability": 0.8},
        designed_sequences=["TESTSEQUENCE"],
        confidence_scores=[0.9],
        warnings=[]
    )
    
    integrator.export_results(result, output_file)
    
    with open(output_file, 'r') as f:
        data = json.load(f)
    
    # Check required fields
    required_fields = [
        'binding_site_score', 'model_used', 'sequence_recovery',
        'structural_plausibility', 'druggability_metrics', 'designed_sequences',
        'confidence_scores', 'warnings', 'timestamp', 'config'
    ]
    
    for field in required_fields:
        assert field in data, f"Missing required field: {field}"
    
    # Check data types
    assert isinstance(data['binding_site_score'], (int, float))
    assert isinstance(data['model_used'], str)
    assert isinstance(data['designed_sequences'], list)
    assert isinstance(data['confidence_scores'], list)
    assert isinstance(data['warnings'], list)
    assert isinstance(data['druggability_metrics'], dict) 