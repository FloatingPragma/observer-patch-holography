"""A finite receiver for the launch phase, with a cap on ordinary separated starts."""
import argparse
import hashlib
from math import exp, log, log1p, sqrt

import numpy as np

from .gateway import HERE, RECEIPT, close, load_json, read, registered_receipt, require, save, settings

PROTOCOL = HERE / "gateway_message_protocol.json"
OUTPUT = HERE / "gateway_message.json"
PROTOCOL_DIGEST = "ee0cffb0138ffb8a7a3bd5b2bbfcf5aa692cc78b160c954a190e9aee06692d14"


def frame(rng, live, cap, admission, click):
    # Conditional on the total successes in cap starts, their positions are
    # a uniform subset. The gaps between them are uniform weak compositions.
    total = int(rng.binomial(cap,admission))
    accepted = min(live,total)
    attempts = cap
    if total >= live:
        gap_fraction = rng.beta(live,total-live+1)
        attempts = live + int(rng.binomial(cap-total,gap_fraction))
    return [attempts,accepted,int(rng.binomial(accepted,click))]


def decode(row, cfg):
    used,admitted,hits = row
    live,cap = cfg["admissions_per_frame"],cfg["maximum_starts"]
    require(all(type(x) is int for x in row) and 0<=hits<=admitted<=live<=cap
            and admitted<=used<=cap and (admitted==live or used==cap), "impossible frame")
    if admitted<live:
        return "erasure"
    return "minus" if hits>=cfg["threshold_minus"] else "plus"


def chernoff(n,cut,p,upper):
    x = cut/n
    require(0<x<1 and 0<p<1 and (x>p if upper else x<p), "uncertified tail direction")
    kl = x*log(x/p) + (1-x)*(log1p(-x)-log1p(-p))
    return max(1e-300,exp(-n*kl))


def replay(record=False):
    digest = hashlib.sha256(PROTOCOL.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    require(digest == PROTOCOL_DIGEST, "registered message protocol differs")
    cfg = load_json(PROTOCOL)
    base = settings()
    live,cap,t = cfg["admissions_per_frame"],cfg["maximum_starts"],cfg["threshold_minus"]
    rng = np.random.default_rng(cfg["seed"])
    receipt = registered_receipt(load_json(RECEIPT))
    rows,bounds = {},{}
    L = log(2*(32+31)*81/base["alpha"])
    for law in ("native","exchange","blocked","dephased"):
        for bit in ("plus","minus"):
            key = f"native_mixed_{cfg['booking']}" if law=="dephased" else f"{law}_{bit}_{cfg['booking']}"
            p = np.array(receipt["weights"][key])/10**13
            admit = sum(p[:-1])
            click = sum(p[:30])/admit
            rows[f"{law}_{bit}"] = [frame(rng,live,cap,admit,click) for _ in range(cfg["frames_per_sign"])]
            if law=="native":
                counts = receipt["counts"][key]
                lo,hi = read(counts,base["shots"],p,L)
                c = 2*L/base["shots"]
                f = sum(counts[:-1])/base["shots"]
                up = f+c/3+c/2+sqrt(c*(f+c/3)+c*c/4)
                admit_lo = f-sqrt(c*up)-c/3-3.2e-12
                erase = chernoff(cap,live-1,admit_lo,False)
                wrong = chernoff(live,t,hi[-1],True) if bit=="plus" else chernoff(live,t-1,lo[-1],False)
                require(wrong+erase < 0.01, "native communication gate failed")
                bounds[bit] = {"wrong_or_erased_upper":float(wrong+erase),"erasure_upper":erase}
    summary = {}
    for name,frames in rows.items():
        bit = name.split("_")[-1]
        decoded = [decode(row,cfg) for row in frames]
        summary[name] = {"wrong":decoded.count("minus" if bit=="plus" else "plus"),
                         "erased":decoded.count("erasure")}
    result = {"frames":rows,"errors":summary,"bounds":bounds}
    if record:
        save(OUTPUT,result)
    else:
        # Even the 1e-300 conservative erasure bound has a meaningful sign.
        close(result, load_json(OUTPUT), absolute=0.0)
    print("Verified all message frames, including competing laws and wrong or erased messages.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", action="store_true")
    replay(parser.parse_args().record)
