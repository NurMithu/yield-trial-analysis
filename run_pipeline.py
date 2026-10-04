"""Run the full analysis:  python run_pipeline.py"""
import sys
sys.path.insert(0, "src")
from yieldtrial.pipeline import run

if __name__ == "__main__":
    res = run()
    print("flow lag:", res["lag"], "s")
    print(res["table"].round(3).to_string(index=False))
    print(res["log"].to_string(index=False))
    print("Moran's I (mixed, +cov):", res["mor"])
    print("variance components:", res["vc"])
    print(res["need"].to_string(index=False))
