import pytest
from selenium import webdriver
from webdriver_manager.chrome import ChromeDriverManager
import os

# Setting up environment variable for Dash to use Chrome
os.environ["DASH_TESTING_BROWSER"] = "chrome"

@pytest.fixture
def driver():
    # Use the WebDriver Manager to automatically manage ChromeDriver
    options = webdriver.ChromeOptions()
    # You can add options like headless, disable GPU, etc.
    options.add_argument("--headless")  # Run Chrome in headless mode
    driver = webdriver.Chrome(ChromeDriverManager().install(), options=options)
    yield driver
    driver.quit()  # Cleanup the driver after each test
