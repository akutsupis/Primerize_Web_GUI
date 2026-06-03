import streamlit as st
from starlette.routing import Route
from starlette.responses import PlainTextResponse

async def robots_txt(request):
    with open("robots.txt", "r") as f:
        content = f.read()
    return PlainTextResponse(content)

# The Starlette app wraps the existing app.py GUI
app = st.App("app.py", routes=[Route("/robots.txt", robots_txt)])
