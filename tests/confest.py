import pytest
from selenium import webdriver
from webdriver_manager.chrome import ChromeDriverManager
import os

os.environ["DASH_TESTING_BROWSER"] = "chrome"

@pytest.fixture
def driver():
    options = webdriver.ChromeOptions()
    options.add_argument("--headless") 
    driver = webdriver.Chrome(ChromeDriverManager().install(), options=options)
    yield driver
    driver.quit() 
